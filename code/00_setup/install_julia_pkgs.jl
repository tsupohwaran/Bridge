#==================================================#
# Install required Julia packages for the project
#==================================================#

using Pkg
ENV["PATH"] = "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:" * ENV["PATH"]

Pkg.activate(pwd())
Pkg.instantiate()

