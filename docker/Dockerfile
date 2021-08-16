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
RUN pip install -q --upgrade pip && pip install -q numpy pandas matplotlib seaborn scikit-hep jupyter && \
    jupyter notebook --generate-config && \
    	echo "c.NotebookApp.ip = '*'" >> ~/.jupyter/jupyter_notebook_config.py && \
	    echo "c.NotebookApp.password = u'sha1:a3e3bea41cca:f13c750230a855de55672091920bd5385ca1ed33'" >> ~/.jupyter/jupyter_notebook_config.py && \
      echo "c.NotebookApp.open_browser = False" >> ~/.jupyter/jupyter_notebook_config.py && \
	    echo "c.NotebookApp.allow_remote_access = True" >> ~/.jupyter/jupyter_notebook_config.py && \
	    echo "c.NotebookApp.port = 8765" >> ~/.jupyter/jupyter_notebook_config.py
ENV PKG_DIR /home/$USER/packages
COPY yadisk.py yadisk.py
COPY install.sh install.sh
RUN sh install.sh
ENV PATH /home/$USER/packages/ISRSolver/bin:/home/$USER/packages/root/bin:$PATH
ENV LD_LIBRARY_PATH /home/$USER/packages/ISRSolver/lib:/home/$USER/packages/root/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH /home/$USER/packages/ISRSolver/lib/python:/home/$USER/packages/root/lib:$PYTHONPATH
ENV ROOTSYS /home/$USER/packages/root
WORKDIR /home/$USER/workdir
EXPOSE 8765

CMD [ "/bin/bash" ]
ENTRYPOINT jupyter notebook --no-browser  --port 8765