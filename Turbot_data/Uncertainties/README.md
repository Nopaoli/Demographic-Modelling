Below are the formula for propagation of uncertainties for a real time estimate from \textit{moments} or \textit{dadi}. You will need the Variance and co variances for all parameter estimates, which you can obtain from the inverse of the Godambe Information Matrix (see python scripts in this folder). \bigskip 

$ \sigma_{T_{real}} = G_t* \mu *\sqrt{\theta^2\sigma^2_{T} + T^2\sigma^2_{\theta}} $ \bigskip

where: \par
$G_t = $ generation time (in this case 3.5 years) \par
$\mu =$  mutation rate for the entire sequence length (in this case, for the "linked" SFS, this is 2660910 bases * 1e-8, assuming a standard mutation rate of $1 \times 10^{-8} $
\end{document}
