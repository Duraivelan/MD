
//  Calculate the force between particles i and j
	
	rij=particle[i].pos-particle[j].pos;
	CorDY=round(rij.comp[1]/box.comp[1]); 
	rij.comp[0]-=shear_rate*box.comp[1]*simu_time*CorDY; 
	rij-=dR;
	rij=rij-shift*xxshift;
	rij2=rij.norm2();

		if (rij2<(r_cut2)) 
		{		
			 if (rij2<(rs2)) {
						    	rij_norm=std::sqrt(rij2);								
												
	// 							exponential potential from PHYSICAL REVIEW E VOLUME 50, NUMBER 3 SEPTEMBER 1994 Browniian dynamics simulations of self-difFusion and shear viscosity of near-hard-sphere colloids
	//							F = n*epsilon*pow(sigma*sigma/rs,n)/rs; //  n is the exponent of the potential function
            					Fij=rij*(Fs/rij_norm);
								particle[i].frc+=Fij;
								particle[j].frc-=Fij;

							//	if(cell[neighborList[i][j][k][m][0]][neighborList[i][j][k][m][1]][neighborList[i][j][k][m][2]][lC2] > cell[i][j][k][lC1]) {
								*p_energy+= phis - phicut + Fs*(rs-rij_norm);
							} 
			else {				
									
				rij2inv=1/rij2;
        		rij6inv 	= 	rij2inv*rij2inv*rij2inv;
        		rij12inv 	= 	rij6inv*rij6inv;
	//			simple potential
	//			F=2*epsilon*(sigma-r)*rx/r;
	//			*p_energy+=2*epsilon*(sigma-r)*(sigma-r);
	
	//			exponential potential from PHYSICAL REVIEW E VOLUME 50, NUMBER 3 SEPTEMBER 1994 Browniian dynamics simulations of self-difFusion and shear viscosity of near-hard-sphere colloids
	//			F = n*epsilon*pow(sigma*sigma/r2,n/2)/r2; //  n is the exponent of the potential function
	//			F = 4*epsilon*(12*pow(sigma*sigma/r2,6)-6*pow(sigma*sigma/r2,3))/r2; 	WCA potential
				F = 4*epsilon*(12*sigma12*rij12inv-6*sigma6*rij6inv)*rij2inv;
        		Fij=rij*F;
								particle[i].frc+=Fij;
								particle[j].frc-=Fij;

				*p_energy+=4*epsilon*(sigma12*rij12inv-sigma6*rij6inv) - phicut;

			} 
			

}
