using Distributions

targetmu = 0.004061023
targetlower = 3.661023/1000.0
targetupper = 4.461023/1000.0


bestdist = Inf
for m=-5.7:0.0001:5.0
	for v=0.01:0.0001:1.0
		global bestdist
		dist = LogNormal(m,v)
		mu = mean(dist)
		lower = quantile(dist, 0.025)
		upper = quantile(dist, 0.975)

		distance = (targetmu-mu)^2.0 + (targetlower-lower)^2.0 + (targetupper-upper)^2.0
		if distance < bestdist
			println(m,"\t",v,"\t",mu, "\t",lower,"\t",upper, "\t",distance)
			bestdist = distance
		end
	end
end