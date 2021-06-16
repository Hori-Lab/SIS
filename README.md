# sis

### Clone

```
git clone --recurse-submodules -j8 git://github.com/naotohori/sis.git
```

### Compile

```
mkdir -p ./build && cd ./build
cmake ..
make -j12
cd ..
ln -s build/sis
```
