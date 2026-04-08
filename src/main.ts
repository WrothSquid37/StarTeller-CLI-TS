import * as np from 'numpy-ts';

interface Cord {
    ra: np.NDArray,
    dec: np.NDArray
}

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

    const partC1 = np.sin(theta_rad).multiply(np.cos(dec_rad)).multiply(ra_rad.add(zeta_rad));
    const partC2 = np.cos(theta_rad).multiply(np.sin(dec_rad));
    const C = partC1.add(partC2);

    const dec_new_rad = np.arcsin(np.clip(C, -1.0, 1.0));
    
    const ra_new_rad = np.arctan2(A, B).add(z_rad);
    
    const ra_new_deg = np.rad2deg(ra_new_rad).mod(360.0);
    const dec_new_deg = np.rad2deg(dec_new_rad);

    const convertCord: Cord = {
        ra: ra_new_deg,
        dec: dec_new_deg
    };

    return convertCord;

}

const jdInput1 = np.array([1,1]);
const jdInput2 = np.array([1,1]);
const lon = np.array([1000.0]);


const result: Cord = precess_cords(jdInput1, jdInput2, lon);

const r1 = result.ra.toArray();
const r2 = result.dec.toArray();



console.log(r1);
console.log(r2);
