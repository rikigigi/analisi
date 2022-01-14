FROM ghcr.io/catthehacker/ubuntu:act-latest
RUN 	curl -sL https://cmake.org/files/v3.16/cmake-3.16.2-Linux-x86_64.sh -o cmakeinstall.sh \
	&& chmod +x cmakeinstall.sh \
	&& ./cmakeinstall.sh --prefix=/usr/local --exclude-subdir \
	&& rm cmakeinstall.sh
