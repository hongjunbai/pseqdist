# pseqdist
A parallel pairwise sequence distance calculator


Supported distances: hamming, identity and similarity.

## Install
git clone the repo, then:
``` shell
> cd pseqdist/
> cargo build --release
> cp target/release/pseqdist ~/bin/
```

## Usage

As an example:
``` shell
> pseqdist -m hamming ./tests/data/toy.fa -o ./test.mat
> diff ./test.mat ./tests/output/toy_seqdist.mat
```

For more infomation
``` shell
> pseqdist -h`
```

## Dependencies
See Cargo.toml
