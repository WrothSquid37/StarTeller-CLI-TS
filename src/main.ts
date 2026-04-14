import * as np from 'numpy-ts';

//Interfaces for storing values that would be returned as tuples in python to make it more clear and consise.
interface Cord {
    ra: np.NDArray,
    dec: np.NDArray
}

interface EquDegrees {
    alt: np.NDArray,
    az: np.NDArray
}

//function unix_timestamp_to_julian_date(ts_unix: np.NDArray) {
  //  return np.
//}

function calc_lst(jd_array: np.NDArray, longitude_deg: number): np.NDArray {
    
	const jd_floor = np.floor(jd_array.subtract(0.5)).add(0.5);
	const day_fraction = jd_array.subtract(jd_floor);

	const T = jd_floor.subtract(2451545.0).divide(36525.0);

    const gmst1 = T.multiply(36000.770053608);
    const gmst2 = np.power(T, 2).multiply(0.000387933);
    const gmst3 = np.power(T, 3).divide(38710000.0);
    const gmst_midnight = gmst1.add(100.46061837).add(gmst2).subtract(gmst3);

    const gmst_deg = day_fraction.multiply(360.98564736629).add(gmst_midnight);

    const lst_deg = gmst_deg.add(longitude_deg).mod(360.0);

	return np.deg2rad(lst_deg);
}

function precess_cords(ra_j2000_deg: np.NDArray, dec_j2000_deg: np.NDArray, jd_target: np.NDArray): Cord {
    
    const t = jd_target.subtract(2451545.0).divide(36525.0);
   
    // Accurate to within a few arcseconds. *split up for clarity and peforming pemdas
    // Rewrite for PEMDAS: (t*num) + ((t^2)*num) - ((t^3)*num)  
    const zeta1 = t.multiply(2306.2181);
    const zeta2 = np.power(t, 2).multiply(0.30188);
    const zeta3 = np.power(t, 3).multiply(0.017998);
    const zeta = zeta1.add(zeta2).subtract(zeta3);

    // Rewrite for PEMDAS: (t*mum) + ((t^2)*num) + ((t^3)*num)
    const z1 = t.multiply(2306.2181);
    const z2 = np.power(t, 2).multiply(1.09468);
    const z3 = np.power(t, 3).multiply(0.018203);
    const z = z1.add(z2).add(z3);

    // Rewrite for PEMDAS: (t*num) - ((t^2)*num) - ((t^3)*num)
    const t1 = t.multiply(2004.3109);
    const t2 = np.power(t, 2).multiply(0.42665);
    const t3 = np.power(t, 3).multiply(0.041833);
    const theta = t1.subtract(t2).subtract(t3);

    // Convert to radians 
    const zeta_rad = np.deg2rad(zeta.divide(3600.0));
    const z_rad = np.deg2rad(z.divide(3600.0));
    const theta_rad = np.deg2rad(theta.divide(3600.0));

    const ra_rad = np.deg2rad(ra_j2000_deg);
    const dec_rad = np.deg2rad(dec_j2000_deg);

    //Precceison formulas
    //Rewrite for formulas for PEMDAS:

    //B1: Cos(theta_rad) * Cos(dec_rad) * Cos(ra_rad + zeta_rad)
    //B2: Sin(theta_rad) * Sin(dec_rad)
    //B: B1 - B2

    //C1: Sin(theta_rad) * Cos(dec_rad) * Cos(ra_rad + zeta_rad)
    //C2: Cos(theta_rad) * Sin(dec_rad)
    //C: C1 + C2

    const A = np.cos(dec_rad).multiply(np.sin(ra_rad.add(zeta_rad)));

    const partB1 = np.cos(theta_rad).multiply(np.cos(dec_rad)).multiply(np.cos(ra_rad.add(zeta_rad)));
    const partB2 = np.sin(theta_rad).multiply(np.sin(dec_rad));
    const B = partB1.subtract(partB2);

    const partC1 = np.sin(theta_rad).multiply(np.cos(dec_rad)).multiply(np.cos(ra_rad.add(zeta_rad)));
    const partC2 = np.cos(theta_rad).multiply(np.sin(dec_rad));
    const C = partC1.add(partC2);

    const dec_new_rad = np.arcsin(np.clip(C, -1.0, 1.0));
    
    const ra_new_rad = np.arctan2(A, B).add(z_rad);
    
    // Convert back to degrees

    const ra_new_deg = np.rad2deg(ra_new_rad).mod(360.0);
    const dec_new_deg = np.rad2deg(dec_new_rad);

    // Convert ra and dec deg into a usable interface for easy access

    const convertCord: Cord = {
        ra: ra_new_deg,
        dec: dec_new_deg
    };

    return convertCord;

}



function equatorial_to_horizontal_deg(ra_deg: np.NDArray, dec_deg: np.NDArray, lst_rad: np.NDArray, lat_rad: np.NDArray, return_azimuth: boolean): np.NDArray; 
function equatorial_to_horizontal_deg(ra_deg: np.NDArray, dec_deg: np.NDArray, lst_rad: np.NDArray, lat_rad: np.NDArray, return_azimuth: boolean): EquDegrees; 
function equatorial_to_horizontal_deg(ra_deg: np.NDArray, dec_deg: np.NDArray, lst_rad: np.NDArray, lat_rad: np.NDArray, return_azimuth: boolean): np.NDArray | EquDegrees {

    /*
    Horizontal altitude (and optionally azimuth) from equatorial coords at a given LST.

    Takes: right ascension, declination, local sidereal time, observer latitude (all radians); return_azimuth boolean
    Returns: Altitude degree and azimuth degree if wanted and only altitude degree if not
    */

    const ra_rad = np.deg2rad(ra_deg);
    const dec_rad = np.deg2rad(dec_deg);

    const ha_rad = lst_rad.subtract(ra_rad);
    
    // SN1: Sin(dec_rad) * Sin(lat_rad)
    // SN2: Cos(dec_rad) * Cos(lat_rad) * Cos(ha_rad)
    // SN: SN1 + SN2
    
    const sn1 = np.sin(dec_rad).multiply(np.sin(lat_rad));
    const sn2 = np.cos(dec_rad).multiply(np.cos(lat_rad)).multiply(np.cos(ha_rad));
    const sin_alt = sn1.add(sn2);


    const alt_rad = np.arcsin(np.clip(sin_alt, -1.0, 1.0));

    const alt_deg = np.rad2deg(alt_rad);

    if (!return_azimuth) {
        return alt_deg;
    }

    const ca = np.cos(alt_rad);
    const cos_alt = np.where(np.abs(ca).less(1e-10), np.array(1e-10), ca);


    const sin_az_neg = np.cos(dec_rad).multiply(-1);

    const sin_az = sin_az_neg.multiply(np.sin(ha_rad)).divide(cos_alt);

    const caz1 = np.sin(dec_rad).subtract(np.sin(lat_rad).multiply(np.sin(alt_rad)));
    const caz2 = np.cos(lat_rad).multiply(cos_alt);
    
    const cos_az = caz1.divide(caz2);

    const az_rad = np.arctan2(sin_az, cos_az);

    const az_deg = np.rad2deg(az_rad).mod(360.0);

    const degs: EquDegrees = {
        alt: alt_deg,
        az: az_deg
    };

    return degs;

}

function sun_equatorial_deg(jd : number) : number;
function sun_equatorial_deg(jd : np.NDArray) : np.NDArray;
function sun_equatorial_deg(jd : np.NDArray | number): number | np.NDArray {
    const is_scalar = np.ndim(jd) == 0;
    if (is_scalar) {
        jd = jd as number;
    }
    
    const n = jd.subtract(2451545.0);
    const g = np.deg2rad(n.multiply(0.9856003).add(357.528).mod(360.0));
    const q = n.multiply(0.985647436).add(280.459).mod(360.0);

    const L1 = np.sin(g).multiply(1.915);
    const L2 = np.sin(g.multiply(2)).multiply(0.020);
    const L = q.add(L1).add(L2);
    
    const L_rad = np.deg2rad(L);

    const e = np.deg2rad(np.subtract(23.439, n.multiply(0.00000036)));
    const ra_rad = np.arctan2(np.cos(e).multiply(np.sin(L_rad)), np.cos(L_rad));
    const dec_rad = np.arcsin(np.sin(e).multiply(np.sin(L_rad)));

    const ra_deg = np.rad2deg(ra_rad).mod(360.0);
    const dec_deg = np.rad2deg(dec_rad);
    if (is_scalar) {
        return 
    }


}
