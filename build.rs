use std::{env, process::Command};

fn main() {
    let out_path_str = env::var("OUT_DIR").unwrap();
    let current_dir = env::current_dir().unwrap().to_str().unwrap().to_string();
    Command::new("sh")
        .arg("-c")
        .arg(&format!("cp -r poa_cpp_code {}", out_path_str))
        .status()
        .unwrap();

    Command::new("sh")
        .arg("-c")
        .arg("/usr/bin/cmake -DCMAKE_BUILD_TYPE:STRING=Release -S./ -B./build -G 'Unix Makefiles'")
        .current_dir(&format!("{}/poa_cpp_code", out_path_str))
        .status()
        .unwrap();

    let _ = Command::new("sh")
        .arg("-c")
        .arg("/usr/bin/cmake --build build/ --config Release --target all -j40 --")
        .current_dir(&format!("{}/poa_cpp_code", out_path_str))
        .status();

    Command::new("sh")
        .arg("-c")
        .arg("/usr/bin/cmake --build build/ --config Release --target all --")
        .current_dir(&format!("{}/poa_cpp_code", out_path_str))
        .status()
        .unwrap();

    let libdir = format!("{}/poa_cpp_code/build/src", out_path_str);
    println!(
        "cargo:rerun-if-changed={}",
        &format!("{}/poa_cpp_code", current_dir)
    );
    println!("cargo:rustc-link-search=native={}", libdir);
    println!("cargo:rustc-link-lib=static=unanimity");
    println!("cargo:rustc-link-lib=static=cc2");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
