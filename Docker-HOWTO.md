Docker HOWTO
============

Basic process:

 1. Build the Docker image
 2. Run the Docker container to execute the desired command

Build the smallk Docker image
-----------------------------

First, make sure you have all submodules and their own submodules:

    git submodule update --init --recursive

Now we can build the image. In the base directory, run this:

    docker build -t smallk .

This will download all dependencies from the Ubuntu repositories, PyPI, GitHub, etc. Everything will be built including
smallk itself. You will end up with a Docker image tagged "smallk".

Running the Docker container
----------------------------

You will need a volume for any input/output data. As an example, you may run the built-in PySmallk tests. The
instructions below assume that your work directory is named `/home/ubuntu`. Replace it with the appropriate name. (The
Docker daemon requires an absolute path for the local volume reference.)

    git clone https://github.com/smallk/smallk_data.git smallk_data
    docker run --volume /home/ubuntu/smallk_data:/data smallk make check PYSMALLK=1 ELEMVER=0.85 DATA_DIR=/data

Here is a breakdown of that Docker command to explain each part:

 - `docker run`: Run a new container from an image
   - `--volume`: Add a volume (persistent storage area) to the container
     - `/home/ubuntu/smallk_data`: Local absolute path that will be exposed within the running container
     - `/data`: Internal path to use within the container
   - `smallk`: Image tag from which to spawn the new container
   - `make check PYSMALLK=1 ELEMVER=0.85`: Run the smallk test suite
     - `DATA_DIR=/data`: Tell the test suite where the local data is stored (from the perspective of the container)
