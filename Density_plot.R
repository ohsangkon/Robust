## SSNSM density line ####
density_ssnsm = function(model, true_epsilon, breaks = 15, freq = F, ylim = NULL,
					main = "", xlab = "Residuals", ylab = "Density"){

	library(sn)	

#	Q = matrix(model$Other_estimates[[component]], ncol = 2)
	Q = model$Other_estimates
	###density function###
	density_Q = 0   # density of NGSM

	for(k in 1 : length(Q)){
		density_Q = density_Q + Q[[k]]$prob * dsn(seq(-20,20,length = 5000), xi = Q[[1]]$xi, omega = Q[[k]]$sigma , alpha = Q[[1]]$lambda)
	}

	par(mar=c(5,5,5,5))
	hist(true_epsilon, breaks = breaks, freq = F, ylim = ylim,
		main = main, xlab = xlab, ylab = ylab, cex.lab = 2.5)
	lines(seq(-20,20,length = 5000), density_Q, col = "red", lwd = 2, lty = 1)
}
