    # R_script.r
# Parse and use command line arguments
# Invoke % Rscript bed_comparison.r $temp_file_name . '_closest_b_to_a.txt'


Args <- commandArgs(TRUE)

my_tab_data <- read.csv(file =Args[1], sep="\t")
my_midpoint_a <- abs((my_tab_data[2]+my_tab_data[3])/2)
my_midpoint_b <- abs((my_tab_data[6]+my_tab_data[7])/2)
my_plot_data <- abs((my_midpoint_a-my_midpoint_b))

#print(sum(my_plot_data[,1]^2), row.names = FALSE)


pdf(paste(Args[2], "_", Args[3], "_", Args[4],"_data_hist.pdf", sep = ""))

list_hist <- hist(my_plot_data[,1], breaks=50, col=Args[5] , xlab=paste("Distance in bp between", "\n", Args[2], "and", Args[3]), main=paste("Histogram of distances between", "\n", Args[2], "and", Args[3]))
text((max(my_plot_data)*0.8),(max(list_hist$counts)*0.9), cex= 1, paste("Sum of the squares = ", "\n", sum(my_plot_data[,1]^2), "\nNumber of locations = ", nrow(my_tab_data), "\nPeak Height = ", Args[4], "\nMean = ", round(mean(my_plot_data[,1]), 0), "\nMedian = ", round(median(my_plot_data[,1]), 0)))


