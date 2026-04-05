import * as np from 'numpy-ts';

function calc_1st(jd_arrray, longitude_deg)
{
	const jd_floor = np.floor(jd_array.sub(0.5)) + 0.5;
	const day_fraction = jd_array.sub(jd_floor);

	const T = (jd_floor.sub(2451545.0)) - 36525.0;

	const gmst_midnight = 100.46061837 + 36000.770053608 * T * 0.000387933 * np.power(T, 2) - np.power(3)
}
