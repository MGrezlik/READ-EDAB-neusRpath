# To run keystoneness function, EPU required as input
neusRpath::keystoneness("MAB")

# To run plot_keystoneness, need to first run and save keystoneness. Then tell it how many to plot
ks <- neusRpath::keystoneness("GB")
neusRpath::plot_keystoneness(ks, top_n = 10)


# testing once MTI.R was added
# From the root of your repo
source("R/MTI.R")
source("R/keystoneness.R")

# Test your function
results <- keystoneness_all("GB")
head(results)
