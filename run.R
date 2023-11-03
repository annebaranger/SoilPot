library(targets)

Sys.setenv("TAR_PROJECT"="data")
tar_visnetwork()
 # tar_make()

# Sys.setenv("TAR_PROJECT"="psi1") # moved in dev, outdated
# tar_visnetwork()
# tar_make()

Sys.setenv("TAR_PROJECT"="psi2")
tar_visnetwork()
# tar_make()

Sys.setenv("TAR_PROJECT"="frost")
tar_visnetwork()
# tar_make()

Sys.setenv("TAR_PROJECT"="traits")
tar_visnetwork()
# tar_make()

Sys.setenv("TAR_PROJECT"="occurence_dataset")
tar_visnetwork()
# tar_make()

Sys.setenv("TAR_PROJECT"="analysis")
tar_visnetwork()
# tar_make()

Sys.setenv("TAR_PROJECT"="modrandom")
tar_visnetwork()
# tar_make()
