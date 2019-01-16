apt-get install curl software-properties-common
curl -sL https://deb.nodesource.com/setup_8.x |bash -
apt-get install nodejs
npm install -g phantomjs-prebuilt --unsafe-perm
echo 'export PATH="${PATH}:/usr/bin/phantomjs"' >> ~/.bashrc
pip install selenium
