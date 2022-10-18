	Lx = 12;
	Ly = 12;
	U = 4;
        N_sites=Lx*Ly;
	N_up=72;
	N_dwn=72;
	dx=0.5;
	dy=sqrt(3.0)/2.0;

	sites=zeros(N_sites,2);
        
	k = 1;
	for i = 1: Ly
		for j = 1: Lx
			if (mod(i,2) == 0)
				x = j-1;
			else
				x = j-1 + dx;
			endif
			y = (1-i)*dy;
			sites(k,1)=x;
			sites(k,2)=y;
			k=k+1;
		end
	end


	H_0 = zeros(2*N_sites,2*N_sites);
	
	for i=1:N_sites
		for j=1:N_sites

			xij = abs(sites(i,1)-sites(j,1));
			xij = min(xij,Lx-xij);
			yij = abs(sites(i,2)-sites(j,2));
			yij = min(yij,Ly*dy-yij);

			r2=xij**2+yij**2;
			if (r2 <= 1.1)
				if (r2 >= 0.9)
					H_0(i,j)=H_0(i,j)-1;
                                        H_0(i+N_sites,j+N_sites)=H_0(i+N_sites,j+N_sites)-1;

				end
			end
		end
	end

   N_part = N_up + N_dwn;

    ueff = min(U,4);
    attempt = 50;
    iteractions = 2000;
    rand('twister',5489);
    emin = 0;
    ibegin = 0;
    for i = 1: attempt
      % create a matrix of random elements
      Hrnd = zeros(2*N_sites,2*N_sites);
        for k=1:2*N_sites
        for l=1:2*N_sites
          Hrnd(k,l)=rand;
        end
      end
      % turn H into a Hermitan matrix
      Hrnd = 0.5*(Hrnd + Hrnd');
      % find the eigenvalues and eigen states of Hrnd
      [evec,evalue] = eig(Hrnd);
      % define the spin up and down sectors of p_trial
      p_trial = horzcat(evec(:,1:N_part));
      % calculate average site densities and energy
        G = p_trial*inv(p_trial'*p_trial)*p_trial';
    
	di_old = zeros(1,N_sites);
	for j=1:N_sites
		di_old(j) = G(j+N_sites,j); 
	end
       e_old=0;
       for j=1:N_part
	       e_old = e_old + evalue(j,j);
       end
      wc = 0;
      it = 1;
      while (wc == 0)
        % construct mean-field Hamiltonians
        HMF = H_0;
	for j=1:N_sites
		%deltaj = ((-1)**j)*ueff*di_old(j);
		deltaj = -ueff*di_old(j);
		HMF(j,j+N_sites) = deltaj;
		HMF(j+N_sites,j) = deltaj;
	end
        % find eigen values and eigen vectors for the spin up sectors
        [evec,evalue]=eig(HMF);
        p_trial = horzcat(evec(:,1:N_part));
        % find eigen values and eigen vectors for the spin down sectors
        % Calculate the Green functions, average densities and e_new
        % calculate the total energy:
        G = p_trial*inv(p_trial'*p_trial)*p_trial';
    
        n_int=0;
	di = zeros(1,N_sites);
	for j=1:N_sites
                di(j) = G(j,j+N_sites);

	end
        e_new = 0;
	for j=1:N_part
		e_new = e_new+evalue(j,j);
	end
	delta_e = abs(e_new-e_old);
        delta_d = 0;
        for k=1:N_sites
	  new = di(k);
	  old = di_old(k);
          delta_d = delta_d + abs(new-old);
        end
        delta_d = delta_d/N_sites;
        if (delta_e <= 10**(-8))
          if (delta_d <= 10**(-4))
	      ibegin = ibegin+1;
              if (ibegin == 1)
                minima(1) = e_new;
                psi_trial = p_trial;
               elseif all (minima - e_new > 10**(-8))
                minima=horzcat(minima,e_new);
                psi_trial = p_trial;
              end
              wc = 1;
          end
        end
	mix = zeros(1,N_sites);
	for k=1:N_sites
        	mix(k) = 0.5*di(k) + 0.5*di_old(k);
        	di_old(k) = di(k);
		di(k) = mix(k);
	end
	e_old=e_new;
	if (it == iteractions)
		wc = 2;
	end
	it = it+1;
      endwhile
      wc
      e_old
      if (wc ==1)
      minima
      end
    end

        psi = psi_trial;
        G = psi*(psi'*psi)*psi';
        n_int=0;
        for i=1:N_sites
                n_int = n_int + G(i,i)*G(i+N_sites,i+N_sites);
                n_int = n_int - G(i,i+N_sites)*G(i+N_sites,i);
        end
        potentialEnergy=n_int*U;
        %% calculate the kinetic energy:
        kineticEnergy=sum(sum(H_0.'.*(G)));
        %% calculate the total energy:
        E_T=potentialEnergy+kineticEnergy;
        E_T
	potentialEnergy
	kineticEnergy

    for i=1:N_up+N_dwn
        for j=1:2*N_sites
                p=psi(j,i);
                save -ascii -append psi.in p
        end
    end

  szsz = zeros(4*Lx + 1,2*Ly + 1);
  nn = zeros (4*Lx + 1,2*Ly + 1);
  count = zeros(4*Lx + 1, 2*Ly + 1);   
  Ns = N_sites;

  for j = 1:N_sites
  for i = 1:N_sites
	  xij=sites(i,1)-sites(j,1);
	  yij=sites(i,2)-sites(j,2);
          %PBC
	  %xij = xij -sidex*round(xij/sidex);
	  %yij = yij -sidey*round(yij/sidey);

          ixij = round(xij/dx) + 2*Lx + 1;
	  iyij = round(yij/dy) + Ly + 1;
	  count(ixij,iyij) = count(ixij,iyij) + 1;
	  if (i==j) 
               c = G(i,i) + G(i+Ns,i+Ns);
                szsz(ixij,iyij) = szsz(ixij,iyij) + c + 2.0*G(i,i+Ns)*G(i+Ns,i) - 2.0*G(i,i)*G(i+Ns,i+Ns);
		nn(ixij,iyij) = nn(ixij,iyij) +  c - c*c;
          else  
                c = G(i,i)*G(j,j)-G(i,j)*G(j,i)+G(i+Ns,i+Ns)*G(j+Ns,j+Ns)-G(i+Ns,j+Ns)*G(j+Ns,i+Ns);
                szsz(ixij,iyij) = szsz(ixij,iyij) + c -G(i,i)*G(j+Ns,j+Ns)+G(i,j+Ns)*G(j+Ns,i)-G(i+Ns,i+Ns)*G(j,j)+G(i+Ns,j)*G(j,i+Ns);
		nn(ixij,iyij) = nn(ixij,iyij) -G(j,i)*G(i,j)-G(j+Ns,i+Ns)*G(i+Ns,j+Ns)-G(j,i+Ns)*G(i+Ns,j)-G(j+Ns,i)*G(i,j+Ns);
          end
  end
  end

f1 = fopen('szsz.dat','w');
for i = 1:4*Lx+1
	for j = 1:2*Ly+1
		if (count(i,j) > 0)
		xij = (i - 2*Lx - 1)*dx;
		yij = (j - Ly - 1)*dy;
    		fprintf(f1,'%d\t%d\t%d\n',xij,yij,szsz(i,j));
		end
	end
end
fclose(f1);

f2 = fopen('nn.dat','w');
for i = 1:4*Lx+1
	for j = 1:2*Ly+1
		if (count(i,j) > 0)
		xij = (i - 2*Lx - 1)*dx;
		yij = (j - Ly - 1)*dy;
    		fprintf(f1,'%d\t%d\t%d\n',xij,yij,nn(i,j));
		end
	end
end
fclose(f2);

bmax = 500;
Ds = zeros(bmax,1);
Dt = zeros(bmax,1);
count = zeros(bmax,1);
dr = 0.1;
sidex = Lx;
sidey = Ly*dy;
for j = 1:Ly-1
	for i = 1+(j-1)*Lx:j*Lx-1
		xij=sites(i+Lx,1)-sites(1+Lx,1);
          	yij=sites(i+Lx,2)-sites(1+Lx,2);

          	%PBC
          	xij = xij -sidex*round(xij/sidex);
          	yij = yij -sidey*round(yij/sidey);
	 	rij = sqrt(xij*xij+yij*yij);
		irij = round(rij/dr)+1;
		s1 = 0.0;
		t1 = 0.0;
		s2 = 0.0;
		t2 = 0.0;
		s3 = 0.0;
		t3 = 0.0;
		% first bound
		i1 = 1;
		i2 = 2;
		j1 = i;
		j2 = i+1;
		s1 = s1+G(i2+Ns,j2+Ns)*G(i1,j1)-G(i2+Ns,j1)*G(i1,j2+Ns);
		s1 = s1+G(i2+Ns,j1+Ns)*G(i1,j2)-G(i2+Ns,j2)*G(i1,j1+Ns);
		s1 = s1+G(i2,j1)*G(i1+Ns,j2+Ns)-G(i2,j2+Ns)*G(i1+Ns,j1);
		s1 = s1+G(i2,j2)*G(i1+Ns,j1+Ns)-G(i2,j1+Ns)*G(i1+Ns,j2);
		t1 = t1+G(i2+Ns,j2+Ns)*G(i1,j1)-G(i2+Ns,j1)*G(i1,j2+Ns);
		t1 = t1-G(i2+Ns,j1+Ns)*G(i1,j2)+G(i2+Ns,j2)*G(i1,j1+Ns);
		t1 = t1-G(i2,j1)*G(i1+Ns,j2+Ns)+G(i2,j2+Ns)*G(i1+Ns,j1);
		t1 = t1-G(i2,j2)*G(i1+Ns,j1+Ns)-G(i2,j1+Ns)*G(i1+Ns,j2);
		% second bound
		i1 = 2;
		i2 = 1+Lx;
		j1 = i+1;
		j2 = i+Lx;
		s2 = s2+G(i2+Ns,j2+Ns)*G(i1,j1)-G(i2+Ns,j1)*G(i1,j2+Ns);
		s2 = s2+G(i2+Ns,j1+Ns)*G(i1,j2)-G(i2+Ns,j2)*G(i1,j1+Ns);
		s2 = s2+G(i2,j1)*G(i1+Ns,j2+Ns)-G(i2,j2+Ns)*G(i1+Ns,j1);
		s2 = s2+G(i2,j2)*G(i1+Ns,j1+Ns)-G(i2,j1+Ns)*G(i1+Ns,j2);
		t2 = t2+G(i2+Ns,j2+Ns)*G(i1,j1)-G(i2+Ns,j1)*G(i1,j2+Ns);
		t2 = t2-G(i2+Ns,j1+Ns)*G(i1,j2)+G(i2+Ns,j2)*G(i1,j1+Ns);
		t2 = t2-G(i2,j1)*G(i1+Ns,j2+Ns)+G(i2,j2+Ns)*G(i1+Ns,j1);
		t2 = t2-G(i2,j2)*G(i1+Ns,j1+Ns)-G(i2,j1+Ns)*G(i1+Ns,j2);
		% thirs bound
		i1 = 1+Lx;
		i2 = 1;
		j1 = i+Lx;
		j2 = i;
		s3 = s3+G(i2+Ns,j2+Ns)*G(i1,j1)-G(i2+Ns,j1)*G(i1,j2+Ns);
		s3 = s3+G(i2+Ns,j1+Ns)*G(i1,j2)-G(i2+Ns,j2)*G(i1,j1+Ns);
		s3 = s3+G(i2,j1)*G(i1+Ns,j2+Ns)-G(i2,j2+Ns)*G(i1+Ns,j1);
		s3 = s3+G(i2,j2)*G(i1+Ns,j1+Ns)-G(i2,j1+Ns)*G(i1+Ns,j2);
		t3 = t3+G(i2+Ns,j2+Ns)*G(i1,j1)-G(i2+Ns,j1)*G(i1,j2+Ns);
		t3 = t3-G(i2+Ns,j1+Ns)*G(i1,j2)+G(i2+Ns,j2)*G(i1,j1+Ns);
		t3 = t3-G(i2,j1)*G(i1+Ns,j2+Ns)+G(i2,j2+Ns)*G(i1+Ns,j1);
		t3 = t3-G(i2,j2)*G(i1+Ns,j1+Ns)-G(i2,j1+Ns)*G(i1+Ns,j2);
	
		Ds(irij) = Ds(irij)+abs(s1)+abs(s2)+abs(s3);
		Dt(irij) = Dt(irij)+abs(t1)+abs(t2)+abs(t3);
		count(irij) = count(irij) + 1.0;
        end
end


f3 = fopen('pair.dat','w');
for i = 1:bmax
    rij = (i-1)*dr;
    if (count(i) > 0)
    	Ds(i) = Ds(i)/(count(i)*6.0);
    	fprintf(f3,'%d\t%d\t%d\n',rij,Ds(i),Dt(i));
    end
end
fclose(f3);

% chiral order parameter
sigmax = [0 1;1 0];
sigmay = [0 -1j;1j 0];
sigmaz = [1 0;0 -1];
% sum over upsidedown triangles

cdt = 0;
k = 0.0;
for j = 1: Ly-1
	for i = (j-1)*Lx + 1  : j*Lx -1
	       for alpha = 1:2
	       for beta = 1:2
	       for gamma = 1:2
	       for delta = 1:2
	       for eps = 1:2
	       for xi = 1:2
		ia = (alpha-1)*Ns + i+Lx;
 		ib = (beta-1)*Ns + i+Lx;
		jg = (gamma-1)*Ns + i;
		jd = (delta-1)*Ns + i;
		ke = (eps-1)*Ns + i+1;
		kx = (xi-1)*Ns + i+1;

		prefac = sigmax(alpha,beta)*(sigmay(gamma,delta)*sigmaz(eps,xi)-sigmaz(gamma,delta)*sigmay(eps,xi));
		prefac = prefac + sigmay(alpha,beta)*(sigmaz(gamma,delta)*sigmax(eps,xi)-sigmax(gamma,delta)*sigmaz(eps,xi));
		prefac = prefac + sigmaz(alpha,beta)*(sigmax(gamma,delta)*sigmay(eps,xi)-sigmay(gamma,delta)*sigmax(eps,xi));

		c = G(ia,ib)*(G(jg,jd)*G(ke,kx)-G(jg,kx)*G(ke,jd));
		c = c + G(ia,jd)*(G(ke,ib)*G(jg,kx)-G(jg,ib)*G(ke,kx));
		c = c + G(ia,kx)*(G(jg,ib)*G(ke,jd)-G(ke,ib)*G(jg,jd));

                cdt = cdt + prefac*c/8;
	       end
	       end
	       end
	       end
	       end
	       end
	       k = k+1.0;
	end
end	
cdt = cdt/k;
cdt

% sum over up triangles

cut = 0;
k = 0.0;
for j = 2: Ly
        for i = (j-1)*Lx + 1  : j*Lx -1
               for alpha = 1:2
               for beta = 1:2
               for gamma = 1:2
               for delta = 1:2
               for eps = 1:2
               for xi = 1:2
                ia = (alpha-1)*Ns + i+1-Lx;
                ib = (beta-1)*Ns + i+1-Lx;
                jg = (gamma-1)*Ns + i+1;
                jd = (delta-1)*Ns + i+1;
                ke = (eps-1)*Ns + i;
                kx = (xi-1)*Ns + i;

                prefac = sigmax(alpha,beta)*(sigmay(gamma,delta)*sigmaz(eps,xi)-sigmaz(gamma,delta)*sigmay(eps,xi));
                prefac = prefac + sigmay(alpha,beta)*(sigmaz(gamma,delta)*sigmax(eps,xi)-sigmax(gamma,delta)*sigmaz(eps,xi));
                prefac = prefac + sigmaz(alpha,beta)*(sigmax(gamma,delta)*sigmay(eps,xi)-sigmay(gamma,delta)*sigmax(eps,xi));

                c = G(ia,ib)*(G(jg,jd)*G(ke,kx)-G(jg,kx)*G(ke,jd));
                c = c + G(ia,jd)*(G(ke,ib)*G(jg,kx)-G(jg,ib)*G(ke,kx)); 
                c = c + G(ia,kx)*(G(jg,ib)*G(ke,jd)-G(ke,ib)*G(jg,jd));  

                cut = cut + prefac*c/8;
               end
               end
               end
               end
               end
               end
               k = k+1.0;
        end
end
cut = cut/k;
cut

chiral = cdt + cut;
chiral

