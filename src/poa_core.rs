use std::ffi::{c_char, c_float, c_int, CStr, CString};

pub trait TSubread {
    fn get_seq(&self) -> &str;
    fn get_cx(&self) -> u8;
}

impl TSubread for &String {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

impl TSubread for String {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

impl TSubread for &str {
    fn get_seq(&self) -> &str {
        self
    }

    fn get_cx(&self) -> u8 {
        3
    }
}

#[repr(C)]
pub struct Subread {
    seq: *const c_char,
    flags: c_int,
}

impl Subread {
    pub fn new(seq: *const c_char, flags: c_int) -> Self {
        Self { seq, flags }
    }
}

#[repr(C)]
pub struct Setting {
    min_identity: c_float,
    match_score: c_int,
    mismatch_score: c_int,
    insertion_score: c_int,
    deletion_score: c_int,
    adjust_strand: c_int,
}

impl Setting {
    pub fn new(
        min_identity: c_float,
        match_score: c_int,
        mismatch_score: c_int,
        insertion_score: c_int,
        deletion_score: c_int,
        adjust_strand: c_int,
    ) -> Self {
        Self {
            min_identity,
            match_score,
            mismatch_score,
            insertion_score,
            deletion_score,
            adjust_strand,
        }
    }

    pub fn set_min_identity(&mut self, min_identity: c_float) -> &mut Self {
        self.min_identity = min_identity;
        self
    }
    pub fn set_match_score(&mut self, match_score: c_int) -> &mut Self {
        self.match_score = match_score;
        self
    }

    pub fn set_mismatch_score(&mut self, mismatch_score: c_int) -> &mut Self {
        self.mismatch_score = mismatch_score;
        self
    }

    pub fn set_insertion_score(&mut self, insertion_score: c_int) -> &mut Self {
        self.insertion_score = insertion_score;
        self
    }

    pub fn set_deletion_score(&mut self, deletion_score: c_int) -> &mut Self {
        self.deletion_score = deletion_score;
        self
    }

    pub fn set_adjust_strand(&mut self, adjust_strand: bool) -> &mut Self {
        if adjust_strand {
            self.adjust_strand = 1;
        } else {
            self.adjust_strand = 0;
        }
        self
    }
}

impl Default for Setting {
    fn default() -> Self {
        Self {
            min_identity: 0.82,
            match_score: 3,
            mismatch_score: -5,
            insertion_score: -2,
            deletion_score: -2,
            adjust_strand: 1,
        }
    }
}

#[repr(C)]
pub struct Result {
    seq: *mut c_char,
    n_passes: usize,
}

impl Result {
    pub fn seq(&self) -> String {
        unsafe { CStr::from_ptr(self.seq).to_str().unwrap().to_string() }
    }

    pub fn n_passes(&self) -> usize {
        self.n_passes
    }
}

impl Drop for Result {
    fn drop(&mut self) {
        unsafe {
            libc::free(self.seq as *mut libc::c_void);
        }
    }
}

extern "C" {
    fn PoaDraftGen(reads: *const Subread, num_sbr: usize, setting: *const Setting) -> Result;
    // fn FreePoaResult(result: Result);
}

pub fn poa_consensus<T: TSubread>(read_infos: &Vec<T>, setting: &Setting) -> (String, usize) {
    let read_infos = read_infos
        .iter()
        .map(|read_info| {
            (
                CString::new(read_info.get_seq()).unwrap(),
                read_info.get_cx() as c_int,
            )
        })
        .collect::<Vec<_>>();
    let reads = read_infos
        .iter()
        .map(|(seq, flag)| Subread::new(seq.as_ptr(), *flag as c_int))
        .collect::<Vec<_>>();

    unsafe {
        let poa_res = PoaDraftGen(reads.as_ptr(), reads.len(), setting);
        // let seq = CString::from_raw(poa_res.seq).to_str().unwrap().to_string();
        let seq = poa_res.seq();
        let npasses = poa_res.n_passes();
        // FreePoaResult(poa_res);
        (seq, npasses)
    }
}
