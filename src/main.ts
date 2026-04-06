import * as np from 'numpy-ts';

function calc_lst(jd_arrray: np.NDArray, longitude_deg: number): np.NDArray
{
	const jd_floor = np.floor(jd_array.subtract(0.5)).add(0.5);
	const day_fraction = jd_array.subtract(jd_floor);

	const T = jd_floor.subtract(2451545.0).divide(36525.0);

	const gmst1 = T.multiply(36000.770053608).multiply(0.000387933).multiply(np.power(T, 2));
	const gmst2 = np.power(T, 3).divide(38710000.0);
	gmst_midnight = 100.46061837 + 36000.770053608 * T + 0.000387933 * T**2 - (T**3) / 38710000.0
	const gmst_deg = gmst_midnight + 360.98564736629 * day_fraction;

	const 1st_deg = gmst_deg + longitude_deg;
	1st_deg = 1st_deg % 360.0;

	return np.deg2rad(1st_deg);
}

