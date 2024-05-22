# TURBULENCE PROJECT

## Introduction

This repository is the result of my master thesis at École Normale Supérieure de Paris, under the supervision of [Alexandros Alexakis](https://www.phys.ens.fr/~alexakis/) and [Emmanuel Dormy](https://www.math.ens.psl.eu/~dormy/), in the field of fluid dynamics, in particular in the study of turbulence. The main goal of the project was to study the dynamics of a 2D fluid flow in a doubly periodic domain, with the presence of a random locallised forcing in the middle of the domain. The project was divided in two main parts: the study of the dynamics of point vortices in a 2D fluid flow by integrating the incompressible Navier-Stokes equations, and the study of the dynamics of a point-vortex model in a 2D fluid flow, in which the _vortex particles_ are driven like an _n_-body system.

## Installation and usage of the code

### Download the codes

To download the codes, clone the repository in your computer:

If you have ssh keys set up in your github account:

```
git clone git@github.com:victorballester7/final-master-thesis.git
cd final-master-thesis
```

If you don't have ssh keys set up in your github account:

```
git clone https://github.com/victorballester7/final-master-thesis.git
cd final-master-thesis
```

### Running the codes

In this repository there are two typse of codes: the ones that can be run in a personal computer and the ones that should be run in a supercomputer.

The two codes that are in the first category are called `pointvortices` and `simple_2DNS`. The two other codes are called `spread_2DNS` and `embarrassed_2DNS`. You can see the details of each code in their respective folders inside the `src` folder.

#### Personal computer codes

The codes are written in Fortran and C++. So you will need the following tools to compile the codes:

- gfortran
- gcc
- make

Additionally, you may want to produce graphics of the results. For that you will need the following libraries:

- gnuplot
- python
- matplotlib
- numpy

First, in order to compile the `simple_2DNS` code, you must install the `fftw-2.1.5` library, located in the `lib` folder. To do so, you can use the following commands:

```
cd lib
tar -xvf fftw-2.1.5.tar
cd fftw-2.1.5
./configure --prefix=/path/where/you/want/to/install/fftw-2.1.5/ --enable-type-prefix
make
make install
```

Then, make sure the respective line where we link the `fftw` library in the `Makefile` of the `simple_2DNS` code is correctly set up. The line should look like this:

```
FFTWLDIR = /usr/local/fftw/2.1.5/lib/
```

To compile and run the code at the same time, you can use the script intended for so. Note that it will also clean the previous files in the folder `data/pointvortices`! So be careful with that, and make a backup of the data if you want to keep it.

```
./pointvortices.sh
```

for the `pointvortices` code, and

```
./simple_2DNS.sh
```

for the `simple_2DNS` code (same thing about the data purging applies here as well).

All the data will be stored in their respective directories inside the `data` folder. If on top of that you want to produce some graphics, you can `gnuplot` and experiment yourself by plotting the data output or use the python scripts to produce some graphics/movies.

For the `pointvortices` code (and again from the root folder of the repository):

```
python src/pointvortices/animate.py
```

or

```
src/pointvortices/animate.sh
```

will produce a movie of the vortices. For the Python version, in order to save the movie in the `videos` folder, just pass any additional argument to the script. For example:

```
python src/pointvortices/animate.py potato
```

For the `simple_2DNS` code (and again from the root folder of the repository) you can do a similar thing:

```
python src/simple_2DNS/movie.py
```

In this case, the movie will automatically be saved under the `videos` folder due to the usually-long time it takes to produce the movie. If you only want to generate the plot slices of the vorticity field, you can use the following command:

```
python src/simple_2DNS/Vis2Db.py
```

and the resulting `.png` files will be stored in the `images` folder.

#### Supercomputer codes

The setup for the supercomputers is a bit more tricky, but it's done in this way in order to, once all finely set up, be faster using the codes. In my project I used two supercomputers: [IDRIS Jean Zay](https://www.idris.fr/jean-zay/) and [MesoPSL](https://wwwmesopsl-new.obspm.fr). The setup for both supercomputers is similar, but they have some small differences.

The aim of this configuration is to be able to work as much as possible locally in your computer with your favourite text editor and configs, and then just use the supercomputer to compile the code and run it.

##### ssh automatic login

In my case I connect to shh from my computer to a computer in the lab, and then from that computer to the supercomputer. To avoid having to type the password every time I connect to the supercomputer, I used the `ssh-copy-id` command. This command copies the public key of your computer to the `authorized_keys` file of the supercomputer.

First create a ssh key in your computer if you don't have one:

```
ssh-keygen -t rsa -b 4096 -C "label_name" -f ~/.ssh/id_rsa
```

Then copy the key to the supercomputer (from the lab computer) using the following command:

```
ssh-copy-id -i ~/.ssh/id_rsa.pub uft62hz@jean-zay.idris.fr
```

where here `uft62hz` is my username in the supercomputer. You can find your username in the supercomputer in the email they sent you when you registered. Since we are using two ssh connections, you will have to do the same for the lab computer:

From your computer:

```
ssh-keygen -t rsa -b 4096 -C "label_name" -f ~/.ssh/id_rsa
ssh-copy-id -i ~/.ssh/id_rsa.pub victor@gauss
```

(`gauss` is the name of my lab computer). Then, in order to access directly to the supercomputer from your computer, instead of doing ssh to your lab computer and then ssh to the supercomputer, you should add the following lines in the `~/.ssh/config` file:

```
Host jean-zay.idris.fr
  ProxyCommand ssh victor@gauss -W %h:%p
```

Note this is the procedure for the `IDRIS Jean Zay` supercomputer. For the `MesoPSL` supercomputer, you first ssh-connect to `vballester@styx.obspm.fr` and from there to `vballester@mesopsl.obspm.fr`. So you will have to adapt the ssh configuration keys accordingly.

##### Installing the fftw library

As in the personal computer case, you will need to install an older version of the `fftw` library. The installation is similar to the one in the personal computer case, but in this case you will have to install the library in the supercomputer. To do so, it's recommended to install it under the `/home/` directory, due to the lack of sudo permissions in the root directory. The installation is done in the same way as in the personal computer case. Don't forget to set the correct path in the `Makefile` of the `spread_2DNS` and `embarrassed_2DNS` codes.

##### IDRIS Jean Zay

Executing the following command in the root folder of the repository

```
./syncIDRIS.sh
```

will sync the code with the supercomputer (in a folder called `CODES` in the `~/` directory). To compile the codes do (change my username `uft62hz` for yours):

```
ssh uft62hz@idris.jean-zay.fr
CODES/yyyymmdd/spread_2DNS/compileIDRIS.sh
```

where `yyyymmdd` is the date of the last syncronization that you did with the supercomputer (same applies for the `embarrassed_2DNS` codes). The code will be compiled in the supercomputer and the necessary files for running the code will be stored in the `$WORK/spread_2DNS/` folder. To run the code, execute:

```
cd $WORK/spread_2DNS/
sbatch jobscriptMPI_IDRIS.slurm
```

And to see that indeed the code is running, run:

```
squeue -u uft62hz
```

##### MesoPSL

The setup for MesoPSL is similar to the one for IDRIS Jean Zay but the paths of the folders differ.

When syncing the code with the supercomputer with

```
./syncMesoPSL.sh
```

the code will be stored in the `~/CODES` folder. To compile the code do (change my username `vballester` for yours):

```
ssh vballester@mesopsl.obspm.fr
CODES/yyyymmdd/spread_2DNS/compileMesoPSL.sh
```

and then to run the code:

```
cd /travail/vballester/spread_2DNS/
sbatch jobscriptMPI_MesoPSL.slurm
```

##### Work even more locally

I created two scripts `openIDRIS.sh` and `openMesoPSL.sh` that mount the important directories of the supercomputer in your computer, in the folders `IDRIS` and `MesoPSL` respectively, located in your desktop directory. You will need though to have the `sshfs` package installed in your computer. To run the scripts, just do:

```
./openIDRIS.sh
./openMesoPSL.sh
```

Note that you will have to change the username of your supercomputer in the scripts. Once you are done with the supercomputer, you can unmount the directories with the following commands:

```
./closeIDRIS.sh
./closeMesoPSL.sh
```

With this you will be able to use your favourite plotter in your computer to plot the results of the supercomputer.

## Some nice results

To give you a glimpse of what you can do with the codes, here are some results that I obtained during my master thesis.

### 2D Navier-Stokes

The following video shows the evolution of the vorticity field of a 2D fluid flow in a doubly periodic domain, for a Reynolds number of 64 and a perturbation region being 8 times smaller than the length of the domain, and the size of the vortices created being 4 times smaller than the size of the perturbation region.

https://github.com/victorballester7/final-master-thesis/assets/78110444/1f183057-01e5-4a4f-bf8f-830e92b06358

### Point vortices
