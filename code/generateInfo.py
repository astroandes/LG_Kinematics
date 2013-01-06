from numpy import *
def generate_info(dmh_filename, pairs_ID_filename, lambda_file, output_file='test.out', hubble=0.73):
    """
    This loads the data from a FOF catalog and a list of pairs IDs. 
    From this information we calculate the distance between the pairs, 
    radial and tangential physical velocities, masses, and angular momentum.

    - "dmh_file" FOF file
    - "pairs_ID_file" Pair indexes file
    - "output_file" File to store the output
    """
    G_GRAV = 4.54E-48 #units of Mpc^3 Msun^-1 s^-2
    KM_TO_MPC = 3.2E-20
    E_UNITS = 1E-36
    L_UNITS = 1.0

    out = open(output_file, 'w')
    outraw = open(output_file+'.raw', 'w')
    halos_prop=loadtxt(dmh_filename) 
    iso_pairs_prop=loadtxt(pairs_ID_filename) 
    lambda_data=loadtxt(lambda_file)

    outraw.write('# M1 M2 x1 y1 z1 x2 y2 z2 v1_x v1_y v1_z v2_x v2_y v2_z\n')
    outraw.write('# Msun Msun Mpc Mpc Mpc Mpc Mpc Mpc km/s km/s kms/ km/s km/s km/s\n')


    #print iso_pairs_prop, iso_pairs_prop.size
    if(iso_pairs_prop.size==2):
        iso_id_1 = empty([1])
        iso_id_2 = empty([1])
        iso_id_1[0] = int_(iso_pairs_prop[0])
        iso_id_2[0] = int_(iso_pairs_prop[1])
    else:
        iso_id_1 = int_(iso_pairs_prop[:,0])
        iso_id_2 = int_(iso_pairs_prop[:,1])

    n_pairs = iso_id_1.size
    print n_pairs
    
    n=0

    for i in range(n_pairs):        
        #pivot on the less massive halo (i.e. the MW)
        id_1 = iso_id_2[i]-1
        id_2 = iso_id_1[i]-1

        props_id_1 = halos_prop[id_1,:]
        props_id_2 = halos_prop[id_2,:]
        
        mass_1 = props_id_1[8]/hubble
        mass_2 = props_id_2[8]/hubble
                
        r_1 = props_id_1[1:4]/hubble
        r_2 = props_id_2[1:4]/hubble

        v_1 = props_id_1[4:7]
        v_2 = props_id_2[4:7]

        v_1_norm = sqrt(sum(v_1*v_1))
        v_2_norm = sqrt(sum(v_2*v_2))

        r_12 = r_2-r_1
        v_12 = v_2-v_1
        
        #eigenvalues
        lambda1_1 = lambda_data[id_1-1,0]
        lambda2_1 = lambda_data[id_1-1,1]
        lambda3_1 = lambda_data[id_1-1,2]

        lambda1_2 = lambda_data[id_2-1,0]
        lambda2_2 = lambda_data[id_2-1,1]
        lambda3_2 = lambda_data[id_2-1,2]

        #relative velocity including Hubble flow
        hubble_flow = hubble*100*r_12
        v_12 = v_12 + hubble_flow 

        #center of mass velocity
        v_cm = (mass_1 * v_1  + mass_2 * v_2)/(mass_1+mass_2)
        v_cm_norm = sqrt(sum(v_cm*v_cm))

        #angular momentum per unit mass
        L_12 = cross(r_12,v_12) 
        L_12_norm = sqrt(sum(L_12 * L_12))

        #energy

        E_kin = (0.5 * sum(v_12 * v_12)) * KM_TO_MPC * KM_TO_MPC # units of Mpc^2 s^-2
        E_pot = - G_GRAV * (mass_1 + mass_2)/ sqrt(sum(r_12*r_12))
        E_kin = E_kin/E_UNITS
        E_pot = E_pot/E_UNITS
        L_12_norm = L_12_norm/L_UNITS

        #separation
        norm_r_12 = sqrt(sum(r_12*r_12))
        total_mass = mass_1 + mass_2
        #print "E", E_kin, E_pot, L_12_norm
        #decompose betwee radial and tangetial velocity
        norm_r_12 = sqrt(sum(r_12*r_12))
        norm_v_12 = sqrt(sum(v_12*v_12))
        unit_r_12 = r_12/norm_r_12

        v_radial_norm = sum(v_12*unit_r_12)
        v_radial = v_radial_norm * unit_r_12

        v_tan = v_12 - v_radial
        v_tan_norm = sqrt(sum(v_tan * v_tan))

        if (v_radial_norm < 0.0 and norm_r_12 < 1.0):
            n = n+1


        if(v_radial_norm<0.0 and norm_r_12 < 1.0):
#            print v_tan_norm, v_radial_norm, norm_v_12,lambda1_1, lambda2_1, lambda3_1, lambda1_2, lambda2_2, lambda3_2
            out.write('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' % 
                      (v_radial_norm, v_tan_norm, v_cm_norm, v_1_norm, 
                       v_2_norm, E_kin, E_pot, L_12_norm, lambda1_1, 
                       lambda2_1, lambda3_1, lambda1_2, lambda2_2, lambda3_2, norm_r_12, total_mass))

            outraw.write('%e %e %f %f %f %f %f %f %f %f %f %f %f %f\n' % 
                      (mass_1, mass_2, r_1[0], r_1[1], r_1[2], r_2[0], r_2[1], r_2[2],
                       v_1[0], v_1[1], v_1[2], v_2[0], v_2[1], v_2[2])) 



        if (int_(props_id_1[0]) != int_(iso_id_2[i])):
            print "Deberian ser iguales", props_id_1[0], iso_id_2[i]
            exit 

    out.close()
    outraw.close()
    print "selected in total", n        

    return
