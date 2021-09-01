FROM archlinux:latest

MAINTAINER Sergei Gribanov <ssgribanov@gmail.com>

RUN pacman -Syy --noconfirm && pacman -S --noconfirm git wget gcc gcc-fortran make cmake nlopt eigen gsl boost \
    binutils libx11 libxpm libxft libxext openssl texlive-most \
    python python-pip vim && \
    pacman -Sc --noconfirm
ENV USER isr
RUN useradd -m -d /home/$USER $USER
USER $USER
WORKDIR /home/$USER
RUN python -m venv isrenv
ENV PATH /home/$USER/isrenv/bin:$PATH
ENV PKG_DIR /home/$USER/packages
COPY docker/yadisk.py yadisk.py
COPY docker/install.sh install.sh
RUN sh install.sh
ENV PATH /home/$USER/packages/ISRSolver/bin:/home/$USER/packages/root/bin:$PATH
ENV LD_LIBRARY_PATH /home/$USER/packages/ISRSolver/lib:/home/$USER/packages/root/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH /home/$USER/packages/ISRSolver/lib/python:/home/$USER/packages/root/lib:$PYTHONPATH
ENV ROOTSYS /home/$USER/packages/root
WORKDIR /home/$USER/workdir
EXPOSE 8765

CMD [ "/bin/bash" ]
ENTRYPOINT jupyter notebook --no-browser  --port 8765