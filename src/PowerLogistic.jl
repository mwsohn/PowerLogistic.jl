module PowerLogistic

using Distributions, Roots, Plots

export powerLogistic, mdpr, mdor, plotmdor

# power for simple logistic regression based on
# http://personal.health.usf.edu/ywu/logistic.pdf
# Hsieh et al. STATISTICS IN MEDICINE
# Statist. Med. 17, 1623-1634 (1998)
#

# created on Sept. 24, 2012

# OR = exp(β) odds ratio
# p1 = Pr(Y = 1)
function SSizeLogisticCon(p1, OR; alpha=0.05, power=0.8)
   β = log(OR)
   za = quantile(Normal(), 1 - alpha/2)
   zb = quantile(Normal(), power)
   return ceil(Int64,(za+zb)^2 / (p1*(1-p1)*β^2))
end

# OR = exp(β) odds ratio
# p1 = Pr(Y = 1)
function powerLogisticCon(n, p1, OR; alpha=0.05)
   β = log(OR)
   za = quantile(Normal(),1-alpha/2)
   return ccdf(Normal(),za-sqrt(n*β^2*p1*(1-p1)))
end

###########################
# logistic regression logit(p) = beta0+ beta1*X
# B=pr(X=1)
# p1 = Pr(D|X=0) # event rate at X=0
# p2 = Pr(D|X=1) # event rate at X=1
# alpha - type I error rate
function SSizeLogisticBin(p1, p2, B; alpha=0.05, power=0.8)
   za = quantile(Normal(),1-alpha/2)
   zb = quantile(Normal(),power)

   p=(1-B)*p1+B*p2
   part1 = za*sqrt(p*(1-p)/B)
   part2 = zb*sqrt( p1*(1-p1)+p2*(1-p2)*(1-B)/B)
   part3 = (p1-p2)^2*(1-B)
   return ceil(Int64,(part1+part2)^2/part3)
end

###########################
# logistic regression logit(p) = beta0+ beta1*X
# B=pr(X=1)
# p1 = Pr(D|X=0) # event rate at X=0
# p2 = Pr(D|X=1) # event rate at X=1
# alpha - type I error rate
function powerLogisticBin(n, p1, p2, B; alpha=0.05)
   za = quantile(Normal(),1-alpha/2)

   p=(1-B)*p1+B*p2
   a = za*sqrt(p*(1-p)/B)
   b = sqrt( p1*(1-p1)+p2*(1-p2)*(1-B)/B)
   myc = (p1-p2)^2*(1-B)
   return cdf(Normal(), (sqrt(n*myc) - a)/b )
end

#############################################################
# public functions
#############################################################

"""
   powerLogistic(n = 0, p1 = .0, p2 = .0, B = .0, OR = .0, alpha = 0.05, power = 0.8)

Estimates power or sample size at type I error = `alpha` and power ≥ `power`.
All options are keyword options. If `n` is not specified (i.e., n = 0), it produces sample size.
Otherwise, it produces power of the sample. Other options are:

- `p1` = Pr(Y = 1) for a continous variable. `p1` = Pr(D|X=0), event rate at X = 0.
- `p2` = Pr(D|X=1), event rate at X = 1.
- `B` = Pr(X=1)
- `alpha`: Type I error (default = 0.05)
- `power`: Power (default = 0.8)

### Example
```jldoctest
julia> powerLogistic(p1=.5,OR=exp(0.405),power=0.95)
317

julia> powerLogistic(n = 317, p1=.5,OR=exp(0.405))
0.950...

julia> powerLogistic(p1 = .4, p2 = .5, B = .5, OR=exp(0.405), alpha = 0.05, power = 0.95)
1281

julia> powerLogistic(n = 1281, p1 = .4, p2 = .5, B = .5, OR=exp(0.405), alpha = 0.05)
0.950...
```
"""
function powerLogistic(;
   n = 0,
   p1 = 0.0,
   p2 = 0.0,
   B = 0.0,
   OR = 0.0,
   alpha = 0.05,
   power = 0.8)

   # return sample size
   if n == 0
      if B == 0.0 # continuous
         if p1 == 0.0 || OR == 0.0
            error("`p1` and `OR` must be greater than zero.")
         end
         return SSizeLogisticCon(p1,OR,alpha = alpha, power = power)
      else
         if p1 == p2 && p1 > 0.0 && p2 > 0.0
            error("`p1` and `p2` must be different and must be greater than zero.")
         end
         return SSizeLogisticBin(p1,p2,B,alpha = alpha, power = power)
      end
   # return power
   else
      if B == 0.0 # continuous
         if p1 == 0.0 || OR == 0.0
            error("`p1` and `OR` must be greater than zero.")
         end
         return powerLogisticCon(n, p1, OR,alpha = alpha)
      else
         if p1 == p2 && p1 > 0.0 && p2 > 0.0
            error("`p1` and `p2` must be different and must be greater than zero.")
         end
         return powerLogisticBin(n, p1, p2, B,alpha = alpha)
      end
   end
end

function mdpr(;n = 0, p1 = 0.0, B = 0.0, alpha = 0.05, power = 0.8)

   if n == 0
      error("`n` is required. `n` is the total sample size.")
   end

   if p1 == 0.0
      error("`p1` is required. `p1` is the event rate at X = 0.")
   end

   if B == 0.0
      error("`B` is required. `B` is the proportion of the sample with X = 1.")
   end

   # find p2 that achieves power greater than the power specified
   return fzero(x->powerLogisticBin(n,p1,x,B,alpha=alpha) - power,p1+1e-10,1.0-1e-10)
end

function mdor(;n = 0, p1 = 0.0, B = 0.0, alpha = 0.05, power = 0.8)

   if n == 0
      error("`n` is required. `n` is the total sample size.")
   end

   if p1 == 0.0
      error("`p1` is required. `p1` is the event rate at X = 0.")
   end

   if B == 0.0
      error("`B` is required. `B` is the proportion of the sample with X = 1.")
   end

   # find p2 that achieves power greater than the power specified
   p2 = fzero(x->powerLogisticBin(n,p1,x,B,alpha=alpha) - power,p1 + 1e-9,1.0-1e-9)
   return oddsratio(p2,p1)
end

function oddsratio(p1,p2)
   return (p1/(1-p1))/(p2/(1-p2))
end

function plotmdor(; nvec::Vector = [], p1::Union{Float64,Vector} = 0.0, B::Union{Float64,Vector} = 0.0, alpha = 0.05, power = 0.8)

   len = length(nvec)
   if len < 2
      error("`nvec` must have 2 or more values for sample size")
   end

   orvec = Array{Float64,2}(length(p1),len)

   for (i,p) in enumerate(p1)
      orvec[i,:] = [mdor(n = x, p1 = p, B = B, alpha = alpha, power = power) for x in nvec]
   end

   # plotly
   plotly()

   xlim_lower = nvec[len] / 40
   orvec = orvec'

   # plot
   nlines = length(p1)

   # plot
   if nlines == 1
      plot(nvec,orvec[:,1],
         marker=(2,.5,:circle,:blue),
         left_margin = 8mm,
         title = "Minimum Detectable Odds Ratio",
         ylabel = "Odds Ratio",
         xlabel = "Sample Size",
         ylims=(1-0.1,maximum(orvec)+0.1),
         xlims=(-xlim_lower,nvec[len]+xlim_lower),
         label=false,
         legend = false
      )
   else
      plot(nvec,orvec[:,1],label = string("p1 = ",p1[1]),marker=(2,.5,:circle))
      for i=2:(nlines-1)
         plot!(nvec,orvec[:,i],label = string("p1 = ",p1[i]),marker=(2,.5,:circle))
      end
      plot!(nvec,orvec[:,nlines],
            label = string("p1 = ",p1[nlines]),
            marker=(2,.5,:circle),
            left_margin = 8mm,
            title = "Minimum Detectable Odds Ratio",
            ylabel = "Odds Ratio",
            xlabel = "Sample Size",
            ylims=(1-0.1,maximum(orvec)+0.1),
            xlims=(-xlim_lower,nvec[len]+xlim_lower)
         )
    end
end

end
