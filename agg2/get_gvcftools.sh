if [ ! -e gvcftools-0.16 ]; then
    wget https://github.com/ctsa/gvcftools/archive/v0.16.tar.gz --no-check-certificate -O gvcftools-0.16.tar.gz
fi
tar -xzf gvcftools-0.16.tar.gz
rm -rf gvcftools-0.16.tar.gz
