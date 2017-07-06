#if !defined(__PHASEUTILS_H)
#define __PHASEUTILS_H

/****************************************************************************/

FINLINE double vmfNorm(const double &kappa) {
//	return m_kapa/(2.0*M_PI*(std::exp(m_kapa) - std::exp(-m_kapa)));
	/* Depending on how sinh is implemented, it may be better to switch for
	 * a numerically more stable implementation, e.g. the one proposed at
	 *  http://www.plunk.org/~hatch/rightway.php
	 * Also, it may be better to calculate the inverse of norm and replace
	 * multiplication in vmfPdf with division.
	 */
	return kappa / (4.0 * M_PI * std::sinh(kappa));
}

/****************************************************************************/

FINLINE double hgPdf(const double &cosTheta, const double &g) {
	return 1 / (4 * M_PI) * (1 - g * g) \
			/ std::pow(1.0 + g * g - 2.0 * g * cosTheta, 1.5);
}

FINLINE double vmfPdf(const double &cosTheta, const double &kappa, \
		const double &norm) {
	return std::exp(kappa * cosTheta) * norm;
}

FINLINE double vmfPdf(const double &cosTheta, const double &kappa) {
	const double norm = kappa / (2.0 * M_PI * (std::exp(kappa) \
			- std::exp(- kappa)));
	return std::exp(kappa * cosTheta) * norm;
}

FINLINE double vmfBackPdf(const double &cosTheta, const double &kappa, \
		const double &norm) {
	return std::exp(- kappa * cosTheta) * norm;
}

FINLINE double vmfBackPdf(const double &cosTheta, const double &kappa) {
	const double norm = kappa / (2.0 * M_PI * (std::exp(kappa) \
			- std::exp(- kappa)));
	return std::exp(- kappa * cosTheta) * norm;
}

/****************************************************************************/

FINLINE double hgInverseCdf(double sample, double g) {
	double sqrTerm = (1 - g * g) / (1 - g + 2 * g * sample);
	return (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
}

FINLINE double vmfInverseCdf(double sample, double kappa) {
	return std::log(sample * std::exp(kappa) + (1.0 - sample) \
			* std::exp(- kappa)) / kappa;
}

FINLINE double vmfBackInverseCdf(double sample, double kappa) {
	return - std::log(sample * std::exp(- kappa) + (1.0 - sample) \
			* std::exp(kappa)) / kappa;
}

#endif /* __PHASEUTILS_H */
