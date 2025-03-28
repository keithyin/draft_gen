use std::{
    collections::HashSet,
    ops::{Deref, DerefMut},
    thread,
    time::Instant,
};

use clap::{self, Parser};
use crossbeam::channel::{Receiver, Sender};
use draft_gen::{poa_consensus, Setting, TSubread};
use gskits::{
    gsbam::bam_record_ext::BamRecordExt,
    pbar::{get_spin_pb, DEFAULT_INTERVAL},
    utils::ScopedTimer,
};
use num_cpus;
use rsedlib::{edlib_align, utils::reverse_complement};
use rust_htslib::bam::{self, header::HeaderRecord, Header, Read, Reader, Record, Records, Writer};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(short = 'i')]
    pub sbr_bam: String,

    #[arg(short = 'n')]
    first_n_channels: Option<usize>,

    #[arg(long = "ed-unify-strand")]
    ed_unify_strand: bool,

    #[arg(long = "minPasses", default_value_t = 3)]
    min_passes: usize,
}

#[derive(Debug, Clone, Default)]
pub struct Subread {
    pub seq: String,
    pub name: String,
    pub channel: u32,
    pub cx: u8,
}

impl Subread {
    pub fn new(seq: String, name: String, channel: u32, cx: u8) -> Self {
        Self {
            seq,
            name,
            channel,
            cx,
        }
    }
}

impl TSubread for Subread {
    fn get_cx(&self) -> u8 {
        self.cx
    }
    fn get_seq(&self) -> &str {
        &self.seq
    }
}

#[derive(Debug, Clone, Default)]
pub struct ChannelSubreads(pub Vec<Subread>);

impl ChannelSubreads {
    pub fn get_best_sbr(&self) -> Option<&str> {
        self.iter()
            .max_by_key(|sbr| {
                let base = if sbr.cx == 3 { 1 } else { 0 };
                base * 10_000_000 + sbr.seq.len()
            })
            .map(|v| v.seq.as_str())
    }
    pub fn num_passes(&self) -> usize {
        self.iter().filter(|sbr| sbr.cx == 3).count()
    }
}

impl Deref for ChannelSubreads {
    type Target = Vec<Subread>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl DerefMut for ChannelSubreads {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<Subread>> for ChannelSubreads {
    fn from(value: Vec<Subread>) -> Self {
        Self(value)
    }
}

struct ConsensusResult {
    channel_id: u32,
    consensus_str: String,
    num_passes: usize,
}

impl ConsensusResult {
    pub fn new(channel_id: u32, consensus_str: String, num_passes: usize) -> Self {
        Self {
            channel_id,
            consensus_str,
            num_passes: num_passes,
        }
    }
}

pub struct ChannelSbrIter<'a> {
    records: Records<'a, bam::Reader>,
    channel_subreads: Option<ChannelSubreads>,
    first_n_channels: usize,
    out_channels: usize,
}

impl<'a> ChannelSbrIter<'a> {
    pub fn new(records: Records<'a, bam::Reader>, first_n_channels: Option<usize>) -> Self {
        Self {
            records,
            channel_subreads: Some(ChannelSubreads::default()),
            first_n_channels: first_n_channels.unwrap_or(usize::MAX),
            out_channels: 0,
        }
    }
}

impl<'a> Iterator for ChannelSbrIter<'a> {
    type Item = ChannelSubreads;
    fn next(&mut self) -> Option<Self::Item> {
        if self.out_channels >= self.first_n_channels {
            return None;
        }

        let res = loop {
            let record = self.records.next();
            if record.is_none() {
                break self.channel_subreads.take();
            }
            let record = record.unwrap().unwrap();
            let record_ext = BamRecordExt::new(&record);
            let seq = record_ext.get_seq();
            let name = record_ext.get_qname();
            let cx = record_ext.get_cx().unwrap();
            let channel_idx = record_ext.get_ch().unwrap();
            let subread = Subread::new(seq, name, channel_idx, cx);

            if let Some(channel_subreads) = &self.channel_subreads {
                if !channel_subreads.is_empty()
                    && channel_subreads.first().unwrap().channel != channel_idx
                {
                    let res = self.channel_subreads.take();
                    self.channel_subreads = Some(ChannelSubreads(vec![subread]));
                    break res;
                } else {
                    self.channel_subreads.as_mut().unwrap().push(subread);
                }
            } else {
                unreachable!();
            }
        };

        let distinct_channels = res
            .as_ref()
            .unwrap_or(&ChannelSubreads::default())
            .iter()
            .map(|sbr| sbr.channel)
            .collect::<HashSet<_>>();
        assert!(distinct_channels.len() <= 1);
        self.out_channels += 1;
        return res;
    }
}

fn bam_reader_worker(
    bam_file: &str,
    sender: Sender<ChannelSubreads>,
    first_n_channels: Option<usize>,
) {
    let mut reader = Reader::from_path(bam_file).unwrap();
    reader.set_threads(2).unwrap();

    let mut scoped_timer = ScopedTimer::new();

    let mut _timer = scoped_timer.perform_timing();

    let mut instant = Instant::now();
    // let records = reader.records();
    let channel_subreads_iter = ChannelSbrIter::new(reader.records(), first_n_channels);

    for channel_subreads in channel_subreads_iter {
        _timer.done_with_cnt(1);
        _timer = scoped_timer.perform_timing();

        if instant.elapsed().as_secs() > 60 {
            println!(
                "bam_reader_speed:{:.4}iter/sec",
                _timer.speed(Some(1000_000_000))
            );
            instant = Instant::now();
        }
        match sender.send(channel_subreads) {
            Ok(_) => {}
            Err(err) => {
                eprintln!("{}", err);
                return;
            }
        }
    }
}

fn consensus_worker(
    recv: Receiver<ChannelSubreads>,
    send: Sender<ConsensusResult>,
    cli: &Cli,
    idx: usize,
) {
    let mut scoped_timer = ScopedTimer::new();
    let mut instant = Instant::now();

    // let ab = unsafe { abpoa_init() };
    for channel_subreads in recv {
        if instant.elapsed().as_secs() > 60 && idx == 0 {
            println!(
                "consensus_worker:{:.4}iter/sec",
                scoped_timer.speed(Some(1000_000_000))
            );
            instant = Instant::now();
        }
        let mut _timer = scoped_timer.perform_timing();

        let cns_res = consensus_core(cli, channel_subreads);
        _timer.done_with_cnt(1);
        if let Some(consensus_res) = cns_res {
            match send.send(consensus_res) {
                Ok(_) => {}
                Err(err) => {
                    eprintln!("{}", err);
                    break;
                }
            }
        }
    }
}

fn consensus_core(cli: &Cli, mut channel_subreads: ChannelSubreads) -> Option<ConsensusResult> {
    channel_subreads = filter_and_sort_subreads(channel_subreads);
    if cli.ed_unify_strand {
        unify_strand(&mut channel_subreads);
    }
    let num_passes = channel_subreads.num_passes();
    let mut setting = Setting::default();
    if cli.ed_unify_strand {
        setting.set_adjust_strand(false);
    }

    let (consensus_seq, num_passes) = if num_passes >= cli.min_passes {
        let msa_result = poa_consensus(&channel_subreads, &setting);
        msa_result
    } else {
        if let Some(seq) = channel_subreads.get_best_sbr() {
            (seq.to_string(), num_passes)
        } else {
            return None;
        }
    };

    let consensus_res =
        ConsensusResult::new(channel_subreads[0].channel, consensus_seq, num_passes);
    return Some(consensus_res);
}

fn consensus_writer(recv: Receiver<ConsensusResult>, fname: &str) {
    // 创建 BAM 头部
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");

    header.push_record(&hd);

    // 创建 BAM 写入器
    let mut writer =
        Writer::from_path(fname, &header, bam::Format::Bam).expect("Failed to create BAM writer");

    writer.set_threads(2).unwrap();
    let pbar = get_spin_pb("writing...".to_string(), DEFAULT_INTERVAL);
    let mut scoped_timer = ScopedTimer::new();
    let mut instant = Instant::now();
    for consensus_res in recv {
        if instant.elapsed().as_secs() > 60 {
            println!(
                "consensus_writer:{:.4}iter/sec",
                scoped_timer.speed(Some(1000_000_000))
            );
            instant = Instant::now();
        }

        let _timer = scoped_timer.perform_timing();

        let write_res = consensus_writer_core(&mut writer, consensus_res);

        pbar.inc(1);

        if let Err(e) = write_res {
            eprintln!("Error writing to BAM file: {}", e);
            break;
        }

        _timer.done_with_cnt(1);
    }
    pbar.finish();
}

fn consensus_writer_core(
    bam_writer: &mut Writer,
    consensus_res: ConsensusResult,
) -> anyhow::Result<()> {
    let mut record = Record::new();
    record.set_qname(format!("channel_{}", consensus_res.channel_id).as_bytes());
    record.set(
        format!("channel/{}/poa", consensus_res.channel_id).as_bytes(),
        None,
        consensus_res.consensus_str.as_bytes(),
        &vec![255; consensus_res.consensus_str.len()],
    );
    record
        .push_aux(b"ch", bam::record::Aux::U32(consensus_res.channel_id))
        .unwrap();
    record
        .push_aux(
            b"np",
            bam::record::Aux::U32(consensus_res.num_passes as u32),
        )
        .unwrap();
    let _ = bam_writer.write(&record)?;
    Ok(())
}

fn filter_and_sort_subreads(channel_subreads: ChannelSubreads) -> ChannelSubreads {
    let mut lengths = channel_subreads
        .iter()
        .filter(|sbr| sbr.cx == 3)
        .map(|sbr| sbr.seq.len())
        .collect::<Vec<_>>();

    lengths.sort();

    let max_len = match lengths.len() {
        0 => 1000000,
        1 => lengths[0] * 2,
        n => lengths[n / 2] + lengths[n / 2 - 1],
    };

    let sbureads = channel_subreads.0;
    let mut filtered_subreads = sbureads
        .into_iter()
        .filter(|sbr| sbr.seq.len() < max_len)
        .collect::<Vec<_>>();
    filtered_subreads.sort_by(|left, right| {
        let median = max_len as f32 / 2.0;
        let left_len = left.seq.len() as f32;
        let left_score = (left_len / median).min(median / left_len);

        let right_len = right.seq.len() as f32;
        let right_score = (right_len / median).min(median / right_len);

        right_score.partial_cmp(&left_score).unwrap()
    });

    filtered_subreads.into()
}

fn unify_strand(channel_subreads: &mut ChannelSubreads) {
    for idx in 1..channel_subreads.len() {
        let anchor = &channel_subreads[0];
        let cur_sbr = &channel_subreads[idx];
        let rc = reverse_complement(cur_sbr.seq.as_bytes());
        let edlib_param = rsedlib::param::EdlibAlignParam::default();
        let res1 = edlib_align(cur_sbr.seq.as_bytes(), anchor.seq.as_bytes(), &edlib_param);
        let res2 = edlib_align(&rc, anchor.seq.as_bytes(), &edlib_param);

        match (res1, res2) {
            (Ok(res1), Ok(res2)) => {
                if res1.edit_distance > res2.edit_distance {
                    channel_subreads[idx].seq = unsafe { String::from_utf8_unchecked(rc) };
                }
            }

            (Err(_), Ok(_)) => {
                channel_subreads[idx].seq = unsafe { String::from_utf8_unchecked(rc) };
            }
            _ => {}
        }
    }
}

fn single_thread(cli: &Cli, oup_bam: &str) {
    let inp_bam = &cli.sbr_bam;
    let mut reader = Reader::from_path(inp_bam).unwrap();
    reader.set_threads(2).unwrap();

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.5");
    hd.push_tag(b"SO", "unknown");

    header.push_record(&hd);

    // 创建 BAM 写入器
    let mut writer =
        Writer::from_path(oup_bam, &header, bam::Format::Bam).expect("Failed to create BAM writer");

    writer.set_threads(2).unwrap();

    let channel_sbr_iter = ChannelSbrIter::new(reader.records(), cli.first_n_channels);

    let pbar = get_spin_pb("processing...".to_string(), DEFAULT_INTERVAL);
    // let ab = unsafe { abpoa_init() };

    for channel_sbr in channel_sbr_iter {
        if let Some(v) = consensus_core(cli, channel_sbr) {
            let _ = consensus_writer_core(&mut writer, v);
        }
        pbar.inc(1);
    }
    pbar.finish();
}

fn multi_threads(cli: &Cli, oup_bam: &str) {
    let inp_bam = &cli.sbr_bam;

    let threads = cli.threads.unwrap_or(num_cpus::get_physical() / 2);
    thread::scope(|thread_scope| {
        let (reader_sender, reader_recv) = crossbeam::channel::bounded(1000);
        thread_scope.spawn(move || {
            bam_reader_worker(&inp_bam, reader_sender, cli.first_n_channels);
        });

        // let mut cnt = 0;
        // let instant = Instant::now();
        // for v in reader_recv {
        //     cnt += 1;
        // }
        // println!("cnt:{cnt}, {}", instant.elapsed().as_secs());

        let (consensus_result_sender, consensus_result_recv) = crossbeam::channel::bounded(1000);
        for thread_idx in 0..threads {
            let reader_recv_ = reader_recv.clone();
            let consensus_result_sender_ = consensus_result_sender.clone();
            thread_scope.spawn(move || {
                consensus_worker(reader_recv_, consensus_result_sender_, cli, thread_idx);
            });
        }
        drop(reader_recv);
        drop(consensus_result_sender);

        consensus_writer(consensus_result_recv, oup_bam);
    });
}

fn main() {
    let cli = Cli::parse();
    let inp_bam = &cli.sbr_bam;
    let oup_bam = format!("{}.poa.bam", inp_bam.rsplit_once(".").unwrap().0);
    let oup_bam = &oup_bam;

    if let Some(num_t) = cli.threads {
        if num_t == 1 {
            single_thread(&cli, oup_bam);
        } else {
            multi_threads(&cli, oup_bam);
        }
    } else {
        multi_threads(&cli, oup_bam);
    }
}

#[cfg(test)]
mod test {
    use crate::{filter_and_sort_subreads, Subread};

    #[test]
    fn test_filter_and_sort_subreads() {
        let subreads = vec![
            Subread::new("ACTTTG".to_string(), "name".to_string(), 0, 3),
            Subread::new("ACTTTTTTG".to_string(), "name".to_string(), 0, 3),
            Subread::new(
                "ACTTTGGGGGACTTTGGGGGACTTTGGGGGACTTTGGGGG".to_string(),
                "name".to_string(),
                0,
                3,
            ),
            Subread::new("CCCCAAAGT".to_string(), "name".to_string(), 0, 3),
            Subread::new("CA".to_string(), "name".to_string(), 0, 3),
        ];

        let res = filter_and_sort_subreads(subreads.into());
        println!("{:?}", res);
    }
}
