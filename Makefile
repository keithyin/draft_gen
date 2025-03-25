build:
	cargo build --release

install:
	cargo build --release
	cp target/release/draft_gen /usr/bin

clean: 
	rm -rf target