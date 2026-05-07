#==================================================#
# load required packages
#==================================================#

using Pkg

using LinearAlgebra,
    Statistics,
    StatFiles, # read and write Stata files
    Optim,
    # Package for data manipulation
    JLD2, # save data using julia native format
    ReadStatTables, # read and write Stata
    DelimitedFiles, # read and write delimited files
    XLSX, # read and write Excel files
    CSV, # read and write csv files
    DataFrames, # same as .dta in Stata
    GLM, # for regression analysis

    # Other packages
    AppleAccelerate, # Apple's BLAS !!!Warning: only for MacOS!!!
    BenchmarkTools, # test the speed of the code
    ClipData, # copy data
    Plots,
    StatsPlots

