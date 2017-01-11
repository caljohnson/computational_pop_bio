# Class 1 first script

# This is a script, in it you can write a list of commands to save
# Anything after the # sign is a comment - R doesn't see it
# This is very useful for you to annotate your code so you can remember what it does later
# For example, to save the code for the plot we just made:
rm(list = ls()) # this clears everything
x = seq(0, 1, by=0.1) # x vector
y = exp(x) # y vector
pdf(file="xyplot.pdf") # create file for saving plot as a pdf
plot(x, y, type="l", xlab="x", ylab="y") # plot
dev.off() # shut down current plot (paired with "pdf" command above)
# save as something (e.g., plotxy.R)
# to run, type source("plotxy.R") at the prompt (make sure you're in the same directory where you saved it)
# now you can go back and edit it anytime

# some rules of scripts:
# - comment as you go (also, it's a good idea to have a "preamble" at the top of a script telling you what it's about)
# - define parameters at the top, never use numbers in the code where you can use parameters
# - any time you write a line of code more than once, write it as a function (now for functions...)
