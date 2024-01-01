Where to download: https://github.com/freebayes/freebayes/issues/770

VERSION=1.3.6
wget https://github.com/freebayes/freebayes/releases/download/v${VERSION}/freebayes-${VERSION}-linux-amd64-static.gz
gunzip freebayes-${VERSION}-linux-amd64-static.gz
mv freebayes-${INSTALLVERSION}-linux-amd64-static freebayes
chmod +x freebayes
wget https://raw.githubusercontent.com/freebayes/freebayes/v${VERSION}/scripts/freebayes-parallel -O freebayes-parallel
chmod +x freebayes-parallel

