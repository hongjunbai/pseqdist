# pseqdist
A parallel pairwise sequence distance calculator


Supported distances: hamming, identity and similarity.

## Install
git clone the repo, then:

> cd pseqdist/

> cargo build --release

> mv target/release/pseqdist ~/bin/

## Usage

As an example:
> pseqdist -m hamming ./tests/data/toy.fa -o ./test.mat

> diff ./test.mat ./tests/output/toy_seqdist.mat

For more infomation
> pseqdist -h

## Dependencies
See "Cargo.toml"
