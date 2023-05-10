# import packages
import pandas as pd
import os


# config files 
configfile: "config/config.yaml"

# search phase -- 1
# refinement phase -- 2
# iterations = config["niterations"]["search"] # should not define within config.yaml, but read from 
evoParams = config["evoParams"]
iterations = pd.read_csv(evoParams).iloc[-1]["iteration_end"]
iterations = int(iterations)

repetitions = 1 ## replace to read from file when needed

# constrain wildcards to only match regular expression
wildcard_constraints:
    iter = "[0-9]+",
    rep = "[0-9]+",
    sim = "[0-9]+"

localrules: all, sampleVAl, runEvolutionary,  gatherResults

# target rule
rule all:
    input:
        "rplot.pdf"
# send execution messages ------------------------------------------------------
# onsuccess:
#     print("Workflow finished, no error")
#     shell("mail -r azadeh.hassanpour@uni-goettingen.de -s \"Workflow finished, no error\" azadeh.hassanpour@uni-goettingen.de < {log}")
# onerror:
#     print("An error occurred")
#     shell("mail -r azadeh.hassanpour@uni-goettingen.de -s \"an error occurred\" azadeh.hassanpour@uni-goettingen.de < {log}")

# define input functions ----------------------------------- 

def aggregate_input(wildcards):
    global evoParams
    n = int(wildcards.iter)
    evoParams_pd = pd.read_csv(evoParams)
    for index, row in evoParams_pd.iterrows():
        if int(row["iteration_start"]) <= n and int(row["iteration_end"]) >= n:
            repetitions = int(row["nrep"])
            break
    if n == 1:
        return("")
    elif n == 2:
        file = checkpoints.sampleVAl.get().output[0]
        file = pd.read_csv(file,index_col = 0)
        files = ["results/iteration_1/results_{sim}_{rep}.RData".format(sim=s,rep=r) for s in file.index for r in range(1,int(repetitions) + 1,1)]
        return(files)
    else:
        file = checkpoints.runEvolutionary.get(iter=n-1).output[0]
        file = pd.read_csv(file,index_col = 0)
        files = ["results/iteration_{iter}/results_{sim}_{rep}.RData".format(iter=n-1,sim=s, rep=r) for s in file.index for r in range(1,int(repetitions) + 1,1)]
        return(files)

def aggregate_input_final(wildcards):
    global iterations
    stop_file = "finish.txt" 
    if os.path.exists(stop_file):
        n = int(pd.read_csv(stop_file, header = None)[0][0])
    else:
        n = int(iterations)
    file = checkpoints.runEvolutionary.get(iter = n).output[0]
    file = pd.read_csv(file,index_col = 0)
    files = ["results/iteration_{iter}/results_{sim}_{rep}.RData".format(iter = n, sim=s, rep=r) for s in file.index for r in range(1,int(repetitions) + 1,1)]
    return(files)



def readParams(wildcards,input):
    file = pd.read_csv(input.parameters,index_col = 0)
    params = file.loc[int(wildcards.sim)].tolist()
    return(params)
    



# rules -----------------------------------------------------
# sample initial values
checkpoint sampleVAl:
    """
    sample starting values for evolutionary algorithm
    """
    input:
        sampleScript = config["sampleScript"]
    output:
        "results/parameters/params_1.csv"
    threads: 1
    log: "results/parameters/params_1.csv.log"
    shell:
        """
        (
        env OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
            Rscript --vanilla {input.sampleScript} {output}
        ) &> {log} 
        """

# search phase
# run first batch of breeding programs
rule runSimulations:
    """
    """
    input:
        parameters = "results/parameters/params_{iter}.csv",
        script = ancient(config["simuScript"])
    output:
        "results/iteration_{iter}/results_{sim}_{rep}.RData"
    params:
        params = readParams
    threads: 2
    resources:
        time = lambda wildcards, attempt: min(attempt,2) * int(eval(config["time"])),
        mem_mb = lambda wildcards, attempt: min(attempt,2) * int(eval(config["memory"]))
    log: "results/iteration_{iter}/results_{sim}_{rep}.RData.log"
    benchmark: "results/iteration_{iter}/results_{sim}_{rep}.RData.bench"
    shell:
        """
        (
        env OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
            Rscript --vanilla {input.script} {output} {wildcards.rep} {params.params} 
        ) &> {log}
        """

# run evolutionary algorithm to derive values for second + batch

checkpoint runEvolutionary:
    """
    """
    input:
        files = aggregate_input,
        script = ancient(config["evoScript"])
    output:
        nextpars = "results/parameters/params_{iter}.csv"
    log:  "results/parameters/params_{iter}.csv.log"
    # wildcard_constraints:
    #     iter = "^[^1]$|^[0-9]{2,}$" 
    threads: 1
    shell:
        """
        (
        env OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
            Rscript --vanilla {input.script} {output.nextpars}  {wildcards.iter} {input.files}
        ) &> {log}
        """

# gatehr final results
rule gatherResults:
    """
    run visualization step
    """
    input:
        files = aggregate_input_final,
        script = config["visualizeParameters"]
    output:
        "rplot.pdf" # che3ck for other output
    threads: 5
    shell:
        """
        echo {input.files} > {output}
        # Rscript --vanilla {input.script} {output}
        """        
#env OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
#    Rscript --vanilla {input.script} {output} {input.files}
        

