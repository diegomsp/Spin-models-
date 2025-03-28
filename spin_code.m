function spin_code(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       I N P U T S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global S cs2 cs2_s J Jq Bmag Bvec D plot_spectra Number_of_states_to_plot_in_spectrum Nc_sim Sim_system_noNc Jmin Jmax Wnc nJ nh heights hmin hmax Jdecay J0 cotunneling cotunneling_system_alone Couplings Vmax nV Temp contrast_max Varying_Couplings
% other variables we need
global   SxR SyR SzR S2 Hs SxRs SyRs SzRs S2s Op Ops S2_system Jr hr JwithNc Er spin_r spin_system_r sz_nc_r sz_system_r CouplingsR Bs namefiless print_eigenstates Jvector tip_polarization branch_purity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S E T   U P   I N T E R N A L   V A R I A B L E S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declare_global_variables(varargin{:})
if (cotunneling == 1)
  V = linspace(-Vmax,Vmax,nV);
  kB = 0.086173324; beta = 1/(kB*Temp);
end
[SxR,SyR,SzR,S2,B] = generate_spin_basis(S);
set_up_variables()
[Bvec] = check_Bvec(Bvec,heights,nJ,nh);
[H0] = construct_spin_hamiltonian_given_coefficients(0,SxR,SyR,SzR,S2,cs2,J,Jq,Bmag,D);
if (cs2_s == 0)
else
    H0 = H0 + cs2_s*S2_system;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C A L C U L A T I O N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%% Nc SIM %%%%%%%%%%%%%%%%%%
if (Nc_sim == 1)

    if (Sim_system_noNc == 1)
        [Cs,Es,spins] = diag_H(Hs,S2s); spins =  0.5*(-1+sqrt(1+4*spins));
        %writematrix(round(Es,6),'Spectra_noNc.txt','Delimiter','tab'); %writematrix(round(spins,6),'Spin_noNc.txt','Delimiter','tab');
        data_out = [[zeros(1,size(Bs,2)),Es']; [zeros(1,size(Bs,2)),spins]; Bs,Cs];  % First row is the energies. Second row is <S^2>. After second row: basis elements and eigenvectors
        %writematrix(round(data_out,6),'System_no_Nc.txt','Delimiter','tab');
        % calculate spin correlations
        %[Scorr] = compute_spin_corr(SxRs,SyRs,SzRs,Cs,'NoNc_');
        if (cotunneling_system_alone == 1)
          [d2IdV2,dIdV,IV] = compute_cotunneling(V,beta,Es,Cs,Ops);
          %writematrix(round([V; IV; [0,dIdV]; [0,0,d2IdV2]],6),'dIdV_system.txt','Delimiter','tab');
          [f] = plot_1D_curve(V(2:end),dIdV,'dIdV_system.jpg','Energy(meV)','dIdV');
          [f] = plot_1D_curve(V(3:end),d2IdV2,'d2IdV2_system.jpg','Energy(meV)','d2IdV2');
        end % end if cotunneling == 1
    end %

    if (heights == 0)

        if ( norm(branch_purity)==0 )
        else 
                Proj = zeros(length(branch_purity),length(Jr),4);
                [H] = construct_spin_hamiltonian_given_coefficients(H0,SxR,SyR,SzR,S2,0*cs2,0*JwithNc,0*Jq,0*Bvec(1,:),0*D);
                [C00] = diag_H(H,S2);  s2m =  0.5*(-1+sqrt(1+4*(diag(C00\S2_system*C00))'  ));
                Mgs = 1 + 2*round(s2m(1)); % multiplicity of the GS
                En0m0 = 1:Mgs;  % states for which the Nc is in the GS and the M in the GS
                if ( abs(round(s2m(Mgs+1)) - round(s2m(1))) < 0.1  )
                    % Nc ex is lower than Mol ex
                    En1m0 = (Mgs+1):3*Mgs; % states for which the Nc is in the ES and the M in the GS
                    Mes = 1 + 2*round(s2m(3*Mgs+1)); % multiplicity of the ES
                    En0m1 = (3*Mgs+1):(3*Mgs + Mes); % states for which the Nc is in the GS and the M in the ES
                    En1m1 = (3*Mgs + Mes + 1):(3*Mgs + 3*Mes); % states for which the Nc is in the ES and the M in the ES
                else 
                    % Nc ex is higher than Mol ex
                    Mes = 1 + 2*round(s2m(Mgs+1)); % multiplicity of the ES
                    En0m1 = (Mgs + 1):(Mgs + Mes); % states for which the Nc is in the GS and the M in the ES
                    En1m0 = (Mgs + Mes +1):(3*Mgs + Mes); % states for which the Nc is in the ES and the M in the GS
                    En1m1 = (3*Mgs + Mes + 1):(3*Mgs + 3*Mes); % states for which the Nc is in the ES and the M in the ES
                end %
        end % end if norm(branch_purity)==0

        for iJnc = 1:length(Jr)
            Jnc = Jr(iJnc);
            [H] = construct_spin_hamiltonian_given_coefficients(H0,SxR,SyR,SzR,S2,0*cs2,Jnc*JwithNc,0*Jq,Bvec(iJnc,:),0*D);
            [C,E,spin] = diag_H(H,S2); Er = [Er; E']; spin_r = [spin_r; 0.5*(-1+sqrt(1+4*spin'))]; spin_system_r = [spin_system_r; 0.5*(-1+sqrt(1+4*(diag(C\S2_system*C))'  )) ]; sz_nc_r = [sz_nc_r; (diag(C\SzR(:,:,1)*C))']; sz_system_r = [sz_system_r; (diag(C\(sum(SzR(:,:,2:end),3))*C))'];
            if (cotunneling == 1)
              if (Varying_Couplings > 0)
                Op = zeros(length(S2),length(S2),3);
                for isite = 1:size(SxR,3)
                  Op(:,:,1) = Op(:,:,1) + CouplingsR(iJnc,isite)*SxR(:,:,isite); Op(:,:,2) = Op(:,:,2) + CouplingsR(iJnc,isite)*SyR(:,:,isite); Op(:,:,3) = Op(:,:,3) + CouplingsR(iJnc,isite)*SzR(:,:,isite);
                end % isite
              end % end if (Varying_Couplings)
              [d2IdV2,dIdV,IV,weightsE] = compute_cotunneling(V,beta,E,C,Op); d2IdV2r(:,iJnc) = [0,0,d2IdV2]; dIdVr(:,iJnc) = [0,dIdV];
            end % end if cotunneling == 1  
            
             if ( norm(branch_purity)==0 )
             else 
                 for is = 1:length(branch_purity)
                     for ii = 1:length(En0m0)
                       Proj(is,iJnc,1) = Proj(is,iJnc,1) + (abs(  C00(:,En0m0(ii))'*C(:,branch_purity(is))  )^2);
                     end % ii
                     for ii = 1:length(En1m0)
                       Proj(is,iJnc,2) = Proj(is,iJnc,2) + (abs(  C00(:,En1m0(ii))'*C(:,branch_purity(is))  )^2);
                     end % ii
                    for ii = 1:length(En0m1)
                       Proj(is,iJnc,3) = Proj(is,iJnc,3) + (abs(  C00(:,En0m1(ii))'*C(:,branch_purity(is))  )^2);
                    end % ii
                    for ii = 1:length(En1m1)
                       Proj(is,iJnc,4) = Proj(is,iJnc,4) + (abs(  C00(:,En1m1(ii))'*C(:,branch_purity(is))  )^2);
                    end % ii
                 end % is
             end % if norm(branch_purity)==0

        end % Jr

        if (cotunneling == 1)
          %writematrix(round([dIdVr],6),strcat('dIdV_map_J_',namefiless,'.txt'),'Delimiter','tab'); writematrix(round([d2IdV2r],6),strcat('d2IdV2_map_J_',namefiless,'.txt'),'Delimiter','tab');
          [f] = plot_2D_map(V,Jr,d2IdV2r,contrast_max, strcat('d2IdV2_map_J_',namefiless,'.pdf'));
        end % cotunneling == 1
        %writematrix(round(Er,6),strcat('Spectra_J_',namefiless,'.txt'),'Delimiter','tab'); %writematrix(round(spin_r,6),'Spin_J.txt','Delimiter','tab'); %writematrix(round(spin_r,6),'Spin_system_J.txt','Delimiter','tab'); %writematrix(round(sz_nc_r,6),'sz_nc_J.txt','Delimiter','tab'); %writematrix(round(sz_system_r,6),'sz_system_J.txt','Delimiter','tab');
        % plot spectra
        if (plot_spectra == 1)
          %[f] = plot_spectra_J(Er,          Jr,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,strcat('Spectra_J_',namefiless,'.pdf'));
          [f] = plot_spectra_J(Er - Er(:,1),Jr,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r, strcat('eigenvalues_with_respect_to_GS_J_',namefiless,'.pdf'));
          %[f] = plot_spectra_h(Er - Er(:,1),Jr,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r, strcat('eigenvalues_with_respect_to_GS_J_',namefiless,'.pdf'));
        end % end if plot_spectra
    
        % plot the projections of the states over spin pure states
        if ( norm(branch_purity)==0 )
        else 
            for is = 1:length(branch_purity)
                f=figure("visible","off"); hold on;
                plot(Jr,Proj(is,:,1),'*','DisplayName','NiCp2 GS, Mol GS'); hold on
                plot(Jr,Proj(is,:,2),'*','DisplayName','NiCp2 ES, Mol GS');
                plot(Jr,Proj(is,:,3),'*','DisplayName','NiCp2 GS, Mol ES');
                plot(Jr,Proj(is,:,4),'*','DisplayName','NiCp2 ES, Mol ES');
                xlabel('J (meV)','FontWeight','bold','FontSize',28)
                %ylabel('Projection','FontWeight','bold','FontSize',28)
                set(gca, "linewidth", 4, "fontsize", 28)
                xticks(0:1:max(Jr));
                lgnd = legend('Location','best','NumColumns',1);
                lgnd.FontSize = 6;
                saveas(f,strcat('Projections_',namefiless,'_',num2str(is),'.pdf'));
            end % is
        end % if norm(branch_purity)==0

    else % now heights == 1

        for ih = 1:length(hr)
            Jnc = Jr(ih);
            [H] = construct_spin_hamiltonian_given_coefficients(H0,SxR,SyR,SzR,S2,0*cs2,Jnc*JwithNc,0*Jq,Bvec(ih,:),0*D);
            [C,E,spin] = diag_H(H,S2); Er = [Er; E']; spin_r = [spin_r; 0.5*(-1+sqrt(1+4*spin'))]; spin_system_r = [spin_system_r; 0.5*(-1+sqrt(1+4*(diag(C\S2_system*C))'  )) ]; sz_nc_r = [sz_nc_r; (diag(C\SzR(:,:,1)*C))']; sz_system_r = [sz_system_r; (diag(C\(sum(SzR(:,:,2:end),3))*C))'];
            if (cotunneling == 1)
              if (Varying_Couplings > 0)
                Op = zeros(length(S2),length(S2),3);
                for isite = 1:size(SxR,3)
                  Op(:,:,1) = Op(:,:,1) + CouplingsR(ih,isite)*SxR(:,:,isite); Op(:,:,2) = Op(:,:,2) + CouplingsR(ih,isite)*SyR(:,:,isite); Op(:,:,3) = Op(:,:,3) + CouplingsR(ih,isite)*SzR(:,:,isite);
                end % isite
              end % end if (Varying_Couplings)
              [d2IdV2,dIdV,IV,weightsE] = compute_cotunneling(V,beta,E,C,Op); d2IdV2r(:,ih) = [0,0,d2IdV2]; dIdVr(:,ih) = [0,dIdV];
            end % end if cotunneling == 1
        end % ih
        if (cotunneling == 1)
          %writematrix(round([dIdVr],6),strcat('dIdV_map_h_',namefiless,'.txt'),'Delimiter','tab'); 
          writematrix(round([d2IdV2r],6),strcat('d2IdV2_map_h_',namefiless,'.txt'),'Delimiter','tab');
          [f] = plot_2D_map(V,hr,fliplr(d2IdV2r),contrast_max,strcat('d2IdV2_map_h_',namefiless,'.pdf'));
          %%% plot map + spectrum
          %[f] = plot_2D_map_plus_spectrum(V,hr,fliplr(d2IdV2r),Er,weightsE,contrast_max,strcat('d2IdV2_map_h_andSpectrum_',namefiless,'.pdf'));
        end % cotunneling == 1
        %writematrix(round(Er,6),strcat('Spectra_h_',namefiless,'.txt'),'Delimiter','tab'); %writematrix(round(spin_r,6),'Spin_h.txt','Delimiter','tab'); %writematrix(round(spin_r,6),'Spin_system_h.txt','Delimiter','tab'); %writematrix(round(sz_nc_r,6),'sz_nc_h.txt','Delimiter','tab'); %writematrix(round(sz_system_r,6),'sz_system_h.txt','Delimiter','tab');
        % plot spectra
        if (plot_spectra == 1)
           %[f] = plot_spectra_h(Er,          hr,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,strcat('Spectra_h_',namefiless,'.pdf'));
            [f] = plot_spectra_h(Er - Er(:,1),hr,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,strcat('eigenvalues_with_respect_to_GS_h_',namefiless,'.pdf'));
        end % end if plot_spectra

    end % end if (heights == 0)


else % here Nc_sim == 0
    %%%%  S P E C T R A   A N D  C O R R E L A T I O N S   O F   A   S P I N    M O D E L   ( N O    N C )
    [C,E,spin] = diag_H(H0,S2);
    data_out = [[zeros(1,size(B,2)),E']; [zeros(1,size(B,2)),spin]; B,C];  % First row is the energies. Second row is <S^2>. After second row: basis elements and eigenvectors
    writematrix(round(data_out,6),'System.txt','Delimiter','tab');
    [Scorr] = compute_spin_corr(SxR,SyR,SzR,C,'');
    if (cotunneling == 1)
          [d2IdV2,dIdV,IV] = compute_cotunneling(V,beta,E,C,Op);
          %writematrix(round([V; IV; [0,dIdV]; [0,0,d2IdV2]],6), strcat('dIdV_system_',namefiless,'.txt'),'Delimiter','tab');
          [f] = plot_1D_curve(V(2:end),dIdV, strcat('dIdV_system_',namefiless,'.pdf'),'Energy(meV)','dIdV');
          [f] = plot_1D_curve(V(3:end),d2IdV2, strcat('d2IdV2_system_',namefiless,'.pdf') ,'Energy(meV)','d2IdV2');
    end % end if cotunneling == 1

    if (print_eigenstates == 1)
        %[f] = print_eigenstates_files(C,B);
        for is = 1:size(C,2) % loop over eigenstates
            namefhere = strcat('State_',num2str(is),namefiless,'.txt');
            fileID = fopen(namefhere,'w');
            for irow = 1:size(B,1)
                % create basisel-string:
                if (S(1)==1)    
                    if (B(irow,1)>0)
                      basisel = 'U';
                    end % end if B(irow,1)>0
                    if (B(irow,1)<0)
                       basisel = 'D';
                    end % end if B(irow,1)>0
                    if (B(irow,1)==0)
                        basisel = '0';
                    end % end if B(irow,1)>0
                else
                    if (B(irow,1)>0)
                      basisel = 'u';
                    end % end if B(irow,1)>0
                    if (B(irow,1)<0)
                       basisel = 'd';
                    end % end if B(irow,1)>0
                end % if (S(1)==1)

                for ic = 2:length(B(irow,:))
                    if (S(ic)==1)    
                        if (B(irow,ic)>0)
                            basisel = strcat(basisel,'U');
                        end % end if B(irow,ic)>0
                        if (B(irow,ic)<0)
                            basisel  = strcat(basisel,'D');
                        end % end if B(irow,ic)>0
                        if (B(irow,ic)==0)
                            basisel  = strcat(basisel,'0');
                        end % end if B(irow,ic)>0
                    else
                        if (B(irow,ic)>0)
                            basisel  = strcat(basisel,'u');
                        end % end if B(irow,ic)>0
                        if (B(irow,ic)<0)
                            basisel = strcat(basisel,'d');
                        end % end if B(irow,ic)>0
                    end % if (S(ic)==1)
                end % ic

                fprintf(fileID,'%s %12.4f\n',basisel,C(irow,is));
            end % irow 
            
            fclose(fileID);
        end % is
    end % end if print_eigenstates == 1

end % (Nc_sim == 1)

end % end function master_script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% S U B R O U T I N E S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[C,E,spin] = diag_H(H,S2)
  [C,E]=eig(H); E = diag(E); [E,iE] = sort(E); C = C(:,iE); spin = (diag(C\S2*C))';
end % end function diag_H

function[Scorr] = compute_spin_corr(SxR,SyR,SzR,C,namef)
 % Compute spin corelations
 % input: 'namef' : name for the file to identify it
   Nsites = size(SxR,3);
   for state = size(C,1):-1:1
       name_file = strcat('Spin_corr_',namef,num2str(state),'.txt');
    for j = 1:Nsites
      for k = j:Nsites
       Scorr(j,k) = (C(:,state)')*(  SxR(:,:,j)* SxR(:,:,k) + SyR(:,:,j)* SyR(:,:,k) + SzR(:,:,j)* SzR(:,:,k)  )*(C(:,state))  - ((C(:,state)')*(  SxR(:,:,j)  )*(C(:,state)))*((C(:,state)')*(  SxR(:,:,k)  )*(C(:,state))) - ((C(:,state)')*(  SyR(:,:,j)  )*(C(:,state)))*((C(:,state)')*(  SyR(:,:,k)  )*(C(:,state))) - ((C(:,state)')*(  SzR(:,:,j)  )*(C(:,state)))*((C(:,state)')*(  SzR(:,:,k)  )*(C(:,state)));
       Scorr(k,j) = Scorr(j,k);
      end % k
    end % j
    writematrix(round(Scorr,6),name_file,'Delimiter','tab');
   end % state
end % end function compute_spin_corr

function[SxR,SyR,SzR,S2,B] = generate_spin_basis(S)
% S: vector of spins of each site in the system.

% Total dimension of the problem
M = prod(2*S+1);
Nsites = length(S); % Nsites: number of sites in the system.


S2=zeros(M); SxT=zeros(M); SyT=zeros(M); SzT=zeros(M); H=zeros(M);
for j = 1:Nsites
    [Sx,Sy,Sz,Sp,Sm]=generic_pauli_matrices(S(j));
    if (j==1)
        Sx_j=Sx; Sy_j=Sy; Sz_j=Sz;
    else
        Sx_j=eye(2*S(1)+1); Sy_j=eye(2*S(1)+1); Sz_j=eye(2*S(1)+1);
    end % if (j==1)
    for k = 2:Nsites
        if (k==j)
            matrx=Sx; matry=Sy; matrz=Sz;
        else
            matrx=eye(2*S(k)+1); matry=eye(2*S(k)+1); matrz=eye(2*S(k)+1);
        end % if (k==j)
        Sx_j = kron(Sx_j,matrx);
        Sy_j = kron(Sy_j,matry);
        Sz_j = kron(Sz_j,matrz);
    end %k
    % ...
    SxR(:,:,j)=Sx_j; SyR(:,:,j)=Sy_j; SzR(:,:,j)=Sz_j;
    SxT = SxT + SxR(:,:,j);
    SyT = SyT + SyR(:,:,j);
    SzT = SzT + SzR(:,:,j);
end %j
S2 = SxT^2 + SyT^2 + SzT^2;

%%% generate basis array of the spin basis:
d = 2*S+1;
for i = 1:M
    ib = i;
    iN = mod(ib - 1,d(Nsites)) + 1;
    ib = (ib - iN)/d(Nsites);
    Shere = S(Nsites):-1:-S(Nsites);
    b(Nsites) = Shere(iN);
    for site = Nsites:-1:2
        ijm = mod(ib,d(site-1)) + 1; ib = (ib - ijm + 1)/d(j-1);
        Shere = S(site-1):-1:-S(site-1);
        b(site-1) = Shere(ijm);
    end % site
    B(i,:) = b;
end % i

end % end function generate_spin_operators_in_full_space


function[Sx,Sy,Sz,Sp,Sm]=generic_pauli_matrices(s)

    n = 2*s+1;

    for a = 1:n
        for b = 1:n
            Sx(a,b) = 0.5*((a==b+1) + (a+1==b))*sqrt((s+1)*(a+b-1)-a*b);
            Sy(a,b) = 1i*0.5*((a==b+1) - (a+1==b))*sqrt((s+1)*(a+b-1)-a*b);
            Sz(a,b) = (s+1-a)*(a==b);
        end % b
    end % a

    Sp = Sx + 1i*Sy;
    Sm = Sx - 1i*Sy;

end % generic pauli matrices

function[H] = construct_spin_hamiltonian_given_coefficients(H,SxR,SyR,SzR,S2,cs2,J,Jq,Bmag,D)
% INPUTS:
% H : Hamiltonian to modify
% Bmag: vector of external mahnetic field
% D : array of magnetic anisotropies
% SxR,SyR,SzR,S2 : spin operators in the total space
Nsites = size(SxR,3);
if (H == 0)
    H = zeros(size(S2,1));
end
H = H + cs2*S2;
for site1 = 1:Nsites
    H = H  +  Bmag(1)*SxR(:,:,site1) + Bmag(2)*SyR(:,:,site1) + Bmag(3)*SzR(:,:,site1);
    H = H  +  D(site1,1)*SxR(:,:,site1)^2 + D(site1,2)*SyR(:,:,site1)^2 + D(site1,3)*SzR(:,:,site1)^2;
  for site2 = 1:Nsites
    H = H + J(site1,site2)*( SxR(:,:,site1)*SxR(:,:,site2)+SyR(:,:,site1)*SyR(:,:,site2)+SzR(:,:,site1)*SzR(:,:,site2) );
    H = H + Jq(site1,site2)*( (SxR(:,:,site1)*SxR(:,:,site2)+SyR(:,:,site1)*SyR(:,:,site2)+SzR(:,:,site1)*SzR(:,:,site2))^2 );
  end % k
end % j
end % end function construct_spin_hamiltonian_given_coefficients


function[f] = plot_spectra_J(Er,x,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,name_fig)
f=figure("visible","off"); hold on;
          for state = 1:Number_of_states_to_plot_in_spectrum
            %name_state = strcat('Sz(Nc)=',num2str(round(sz_nc_r(1,state),3)),' S2=',num2str(round(spin_system_r(1,state),3)),' Sz=',num2str(round(sz_system_r(1,state),3)));
            name_state = strcat('Sz(Nc)=',num2str(sz_nc_r(1,state)),' S2=',num2str(spin_system_r(1,state)),' Sz=',num2str(sz_system_r(1,state)));
            plot(x,Er(:,state),'+','LineWidth',4,'DisplayName',name_state);
           end % state

%lgd = legend; lgd.FontSize = 6; lgd.Title.String = 'States'; lgd.NumColumns = 2; lgd.Location = 'northwest';

xlabel('J (meV)','FontWeight','bold','FontSize',28)
%xlabel('height (A)','FontWeight','bold','FontSize',28)
ylabel('Eigenvalues (meV)','FontWeight','bold','FontSize',28)
set(gca, "linewidth", 4, "fontsize", 28)
xticks([0:1:max(x)]);
saveas(f,name_fig);
end % end function plot_spectra_J

function[f] = plot_spectra_h(Er,x,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,name_fig)
f=figure("visible","off"); hold on;

          for state = 1:Number_of_states_to_plot_in_spectrum
            %colormulti='r';
            %if (abs(sval(state,ij)-0.5)<=0.01); colormulti='b'; end %duplet
            %if (abs(sval(state,ij)-1.5)<=0.01); colormulti='g'; end %quadruplet
            %plot(x(ix),E(state,ij)-E(1,ij),'+','LineWidth',4,'Color',colormulti)
            %plot(x,Er(state,:),'+','LineWidth',4,'Color',colormulti)
            name_state = strcat('Sz(Nc)=',num2str(round(sz_nc_r(1,state),3)),' S2=',num2str(round(spin_system_r(1,state),3)),' Sz=',num2str(round(sz_system_r(1,state),3)));
            plot(x,Er(:,state),'+','LineWidth',4,'DisplayName',name_state);
           end % state

%lgd = legend; lgd.FontSize = 6; lgd.Title.String = 'States'; lgd.NumColumns = 2; lgd.Location = 'northwest';

xlabel('height (A)','FontWeight','bold','FontSize',28)
ylabel('Eigenvalues (meV)','FontWeight','bold','FontSize',28)
set(gca, "linewidth", 4, "fontsize", 28)
saveas(f,name_fig);
end % end function plot_spectra

function[d2IdV2,dIdV,IV,weightsE] = compute_cotunneling(V,beta,E,C,Op)
  [populations] = compute_populations(E,beta);
  IV = 0*V(1:end); dIdV = 0*V(2:end); d2IdV2 = 0*V(3:end);
        for n = 1:length(E)
            for m = 1:length(E)
                SW = 0;
                for iop = 1:size(Op,3)
                  SW = SW + abs(C(:,n)'*Op(:,:,iop)*C(:,m))^2;                                     % spectral weight
                end % iop
                weightsE(n,m) = SW;
                Iwp    =   Gfunction((V-(E(m)-E(n))),beta) + Gfunction((V+(E(m)-E(n))),-beta);
                IV     =   IV +  populations(n)*SW*Iwp;
                dIwp   =   diff(Iwp)./diff(V);
                dIdV   =   dIdV +  populations(n)*SW*dIwp;
                d2Iwp  =   diff(dIwp)./diff(V(2:end));
                d2IdV2 =   d2IdV2 +  populations(n)*SW*d2Iwp;
            end % m
        end % n
end % end function compute_cotunneling

function[G] = Gfunction(w, beta)
    G = w./ (1 - exp(-beta*w));
    for l = 1:length(G)
      if (isnan(G(l)) )
        G(l) = 1/beta;
      end % end if (nan(G(l))
    end % l
end %function

function[p] = compute_populations(E,beta)
      Q = 0;
      for n = 1:length(E)
           p(n) = exp(-beta*(E(n)-E(1)));
           Q = Q + p(n);
       end % n
       p = p/Q;
end % end function compute_populations

function[Bvec] = check_Bvec(Bvec,heights,nJ,nh)
if (norm(Bvec)==0)
  if (heights==0)
    Bvec = zeros(nJ,3);
  else
    Bvec = zeros(nh,3);
  end
else % if (norm(Bvec) > 0)
  if (heights==0)
   if (nJ == size(Bvec,1))
   else
     display('Error! Bvec has to have the dimensions of Jr!')
   end
  else
    if (nh == size(Bvec,1))
   else
     display('Error! Bvec has to have the dimensions of hr!')
   end
  end
end % end if (norm(Bvec)==0)
end % end function check_Bvec

function[f] = plot_1D_curve(x,y,name_fig,labelx,labely)
  f=figure("visible","off"); hold on;
  plot(x,y,'LineWidth',4,'DisplayName',name_fig);
  xlabel(labelx,'FontWeight','bold','FontSize',28)
  ylabel(labely,'FontWeight','bold','FontSize',28)
  set(gca, "linewidth", 4, "fontsize", 28)
  saveas(f,name_fig);
end % plot_1D_curve




function[f] = plot_2D_map(x,y,Z,contrast_max,name_map);
   f=figure("visible","off"); hold on;
   imagesc([min(x) max(x)], [min(y) max(y)], flipud(Z'));
   caxis([-contrast_max contrast_max]);
   saveas(f,name_map);
   % Adjust the aspect ratio
   axis xy;
   axis tight;
   pbaspect([1 1 1]); % This keeps the aspect ratio 'auto', similar to matplotlib's 'aspect' parameter
   colorbar;
end % end function plot_2D_map

function[f] = plot_2D_map_NEW(x,y,Z,contrast_max,name_map);
   % Generate or load a 'coolwarm' colormap
n = 256; % Number of colors in the colormap
coolwarm = flipud([linspace(0,1,n)', linspace(1,0,n)', ones(n,1)]); % Approximation of 'coolwarm' (red to blue)

% Alternatively, use colormap from a Python-like external file
% coolwarm = load('coolwarm_colormap.txt'); % Precomputed colormap file

% Flip the colormap to get 'coolwarm_r'
coolwarm_r = flipud(coolwarm);

% Apply to heatmap
f = figure("visible", "off");
hold on;

% Your 2D heatmap
imagesc([min(x) max(x)], [min(y) max(y)], flipud(Z'));
colormap(coolwarm_r); % Use the reversed colormap
caxis([-contrast_max contrast_max]);

% Adjust visualization
colorbar;
axis xy;
axis tight;
pbaspect([1 1 1]); % Aspect ratio 'auto' similar to matplotlib
saveas(f, name_map);
end % end function plot_2D_map

function[f] = plot_2D_map_plus_spectrum(x,y,Z,E,weightsE,contrast_max,name_map)
   f=figure("visible","off"); hold on;
   imagesc([min(x) max(x)], [min(y) max(y)], flipud(Z'));
   caxis([-contrast_max contrast_max]);
   % Adjust the aspect ratio
   axis xy;
   axis tight;
   pbaspect([1 1 1]); % This keeps the aspect ratio 'auto', similar to matplotlib's 'aspect' parameter
   colorbar;
   %%% plot spectral lines
   weightsE = abs(weightsE)/max(max(abs(weightsE(1,:))));
   for i = 2:size(E,2)
       if ( max( E(:,i) - E(:,1) ) < 1.1*max(x)   )
        %colorhere = weightsE(1,i)*[.7  .7  .7];
        plot(  ( E(:,i) - E(:,1) )', y , '--'  , 'LineWidth',4, 'color', 'k' );
       end % end if
   end % i
   set(gca,'YTick',[]);
   set(gca,'YTickLabel',[]);
   box off
   saveas(f,name_map);
end % end function plot_2D_map

function declare_global_variables(varargin)
    % Declare global variables
    % variables to parse (inputs)
    global S cs2 J Jq Bmag Bvec D plot_spectra Number_of_states_to_plot_in_spectrum cs2_s
    global Nc_sim Sim_system_noNc Jmin Jmax Wnc nJ nh hmin hmax Jdecay J0 heights
    global cotunneling cotunneling_system_alone Couplings Vmax nV Temp contrast_max Varying_Couplings
    % other variables we need
    global   SxR SyR SzR S2 Hs SxRs SyRs SzRs S2s Op Ops S2_system
    global Jr hr JwithNc
    global Er spin_r spin_system_r sz_nc_r sz_system_r namefiless print_eigenstates
    global Jvector tip_polarization branch_purity

    % Create an input parser object
    p = inputParser;

    % Define default values
    default_S = [1, 0.5];
    default_cs2 = 0;
    default_cs2_s = 0;
    default_J = [0, 0; 0, 0];
    default_Jq = [0, 0; 0, 0];
    default_Bmag = [0, 0, 0];
    default_Bvec = [0, 0; 0, 0];
    default_D = [0, 0, 0; 0, 0, 0];
    default_plot_spectra = 1;
    default_Number_of_states_to_plot_in_spectrum = 0;
    default_Nc_sim = 0;
    default_Sim_system_noNc = 0;
    default_Jmin = 0.01;
    default_Jmax = 2.00;
    default_nJ = 100;
    default_Wnc = [0,0];
    default_heights = 0;
    default_hmin = 0.01;
    default_hmax = 2.00;
    default_nh = 100;
    default_Jdecay = 0.75;
    default_J0 = 2;
    default_cotunneling = 1;
    default_cotunneling_system_alone = 1;
    default_Couplings = [0,0; 0, 0];
    default_Vmax = 20;
    default_nV = 100;
    default_Temp = 1;
    default_contrast_max = 20;
    default_Varying_Couplings = 0;
    default_namefiless = '';
    default_print_eigenstates = 0;
    default_Jvector = [0,0];
    default_tip_polarization = 0;
    default_branch_purity = [0,0];



    % Define validation functions (if needed)
    validate_S = @(x) isnumeric(x) && isvector(x);
    validate_cs2 = @(x) isnumeric(x) && isscalar(x);
    validate_cs2_s = @(x) isnumeric(x) && isscalar(x);
    validate_J = @(x) isnumeric(x) && ismatrix(x);
    validate_Jq = @(x) isnumeric(x) && ismatrix(x);
    validate_Bmag = @(x) isnumeric(x) && isvector(x);
    validate_Bvec = @(x) isnumeric(x) && ismatrix(x);
    validate_D = @(x) isnumeric(x) && ismatrix(x);
    validate_plot_spectra = @(x) isnumeric(x) && isscalar(x);
    validate_Number_of_states_to_plot_in_spectrum = @(x) isnumeric(x) && isscalar(x);
    validate_Nc_sim = @(x) isnumeric(x) && isscalar(x);
    validate_Sim_system_noNc = @(x) isnumeric(x) && isscalar(x);
    validate_Jmin = @(x) isnumeric(x) && isscalar(x);
    validate_Jmax = @(x) isnumeric(x) && isscalar(x);
    validate_nJ = @(x) isnumeric(x) && isscalar(x);
    validate_Wnc = @(x) isnumeric(x) && isvector(x);
    validate_heights = @(x) isnumeric(x) && isscalar(x);
    validate_hmin = @(x) isnumeric(x) && isscalar(x);
    validate_hmax = @(x) isnumeric(x) && isscalar(x);
    validate_nh = @(x) isnumeric(x) && isscalar(x);
    validate_Jdecay = @(x) isnumeric(x) && isvector(x);
    validate_J0 = @(x) isnumeric(x) && isscalar(x);
    validate_cotunneling = @(x) isnumeric(x) && isscalar(x);
    validate_cotunneling_system_alone = @(x) isnumeric(x) && isscalar(x);
    validate_Couplings = @(x) isnumeric(x) && ismatrix(x);
    validate_Vmax = @(x) isnumeric(x) && isscalar(x);
    validate_nV = @(x) isnumeric(x) && isscalar(x);
    validate_Temp = @(x) isnumeric(x) && isscalar(x);
    validate_contrast_max = @(x) isnumeric(x) && isscalar(x);
    validate_Varying_Couplings = @(x) isnumeric(x) && isscalar(x);
    validate_namefiless = @ischar;
    validate_print_eigenstates = @(x) isnumeric(x) && isscalar(x);
    validate_Jvector = @(x) isnumeric(x) && isvector(x);
    validate_tip_polarization = @(x) isnumeric(x) && isscalar(x);
    validate_branch_purity = @(x) isnumeric(x) && isvector(x);

    % Add optional parameters to the input parser
    addParameter(p, 'Spins', default_S, validate_S);
    addParameter(p, 'cs2', default_cs2, validate_cs2);
    addParameter(p, 'cs2_system', default_cs2_s, validate_cs2_s);
    addParameter(p, 'Jmatrix', default_J, validate_J);
    addParameter(p, 'Jqmatrix', default_Jq, validate_Jq);
    addParameter(p, 'Magnetic_field', default_Bmag, validate_Bmag);
    addParameter(p, 'Varying_Magnetic_field', default_Bvec, validate_Bvec);
    addParameter(p, 'Anisotropies', default_D, validate_D);
    addParameter(p, 'Plot_spectra', default_plot_spectra, validate_plot_spectra);
    addParameter(p, 'Number_of_states_to_plot_in_spectrum', default_Number_of_states_to_plot_in_spectrum, validate_Number_of_states_to_plot_in_spectrum);
    addParameter(p, 'Nc_simulation', default_Nc_sim, validate_Nc_sim);
    addParameter(p, 'Simulation_without_Nc', default_Sim_system_noNc, validate_Sim_system_noNc);
    addParameter(p, 'Jmin', default_Jmin, validate_Jmin);
    addParameter(p, 'Jmax', default_Jmax, validate_Jmax);
    addParameter(p, 'steps_J', default_nJ, validate_nJ);
    addParameter(p, 'Nc_couplings', default_Wnc, validate_Wnc);
    addParameter(p, 'height', default_heights, validate_heights);
    addParameter(p, 'hmin', default_hmin, validate_hmin);
    addParameter(p, 'hmax', default_hmax, validate_hmax);
    addParameter(p, 'steps_h', default_nh, validate_nh);
    addParameter(p, 'Decay', default_Jdecay, validate_Jdecay);
    addParameter(p, 'Maximal_J', default_J0, validate_J0);
    addParameter(p, 'do_cotunneling', default_cotunneling, validate_cotunneling);
    addParameter(p, 'do_cotunneling_system_alone', default_cotunneling_system_alone, validate_cotunneling_system_alone);
    addParameter(p, 'Couplings', default_Couplings, validate_Couplings);
    addParameter(p, 'Vmax', default_Vmax, validate_Vmax);
    addParameter(p, 'steps_V', default_nV, validate_nV);
    addParameter(p, 'Temperature', default_Temp, validate_Temp);
    addParameter(p, 'Contrast', default_contrast_max, validate_contrast_max);
    addParameter(p, 'Varying_Couplings', default_Varying_Couplings, validate_Varying_Couplings);
    addParameter(p, 'namefiless', default_namefiless, validate_namefiless);
    addParameter(p, 'print_eigenstates', default_print_eigenstates, validate_print_eigenstates);
    addParameter(p, 'Jvector', default_Jvector, validate_Jvector);
    addParameter(p, 'tip_polarization', default_tip_polarization, validate_tip_polarization);
    addParameter(p, 'branch_purity', default_branch_purity, validate_branch_purity);



    % Parse the input
    parse(p, varargin{:});

    % Access the values of the optional parameters
    S                                         =          p.Results.Spins;
    cs2                                       =          p.Results.cs2;
    cs2_s                                     =          p.Results.cs2_system;
    J                                         =          p.Results.Jmatrix;
    Jq                                        =          p.Results.Jqmatrix;
    Bmag                                      =          p.Results.Magnetic_field;
    Bvec                                      =          p.Results.Varying_Magnetic_field;
    D                                         =          p.Results.Anisotropies;
    plot_spectra                              =          p.Results.Plot_spectra;
    Number_of_states_to_plot_in_spectrum      =          p.Results.Number_of_states_to_plot_in_spectrum;
    Nc_sim                                    =          p.Results.Nc_simulation;
    Sim_system_noNc                           =          p.Results.Simulation_without_Nc;
    Jmin                                      =          p.Results.Jmin;
    Jmax                                      =          p.Results.Jmax;
    nJ                                        =          p.Results.steps_J;
    Wnc                                       =          p.Results.Nc_couplings;
    heights                                   =          p.Results.height;
    hmin                                      =          p.Results.hmin;
    hmax                                      =          p.Results.hmax;
    nh                                        =          p.Results.steps_h;
    Jdecay                                    =          p.Results.Decay;
    J0                                        =          p.Results.Maximal_J;
    cotunneling                               =          p.Results.do_cotunneling;
    cotunneling_system_alone                  =          p.Results.do_cotunneling_system_alone;
    Couplings                                 =          p.Results.Couplings;
    Vmax                                      =          p.Results.Vmax;
    nV                                        =          p.Results.steps_V;
    Temp                                      =          p.Results.Temperature;
    contrast_max                              =          p.Results.Contrast;
    Varying_Couplings                         =          p.Results.Varying_Couplings;
    namefiless                                =          p.Results.namefiless;
    print_eigenstates                         =          p.Results.print_eigenstates;
    Jvector                                   =          p.Results.Jvector;
    tip_polarization                          =          p.Results.tip_polarization;
    branch_purity                             =          p.Results.branch_purity;


end % end function declare_global_variables


function set_up_variables()
  global S cs2 cs2_s J Jq Bmag Bvec D plot_spectra Number_of_states_to_plot_in_spectrum Nc_sim Sim_system_noNc Jmin Jmax Wnc nJ nh heights hmin hmax Jdecay J0 cotunneling cotunneling_system_alone Couplings Vmax nV Temp contrast_max Varying_Couplings
  % other variables we need
  global   SxR SyR SzR S2 Hs SxRs SyRs SzRs S2s Op Ops S2_system Jr hr JwithNc Er spin_r spin_system_r sz_nc_r sz_system_r CouplingsR Bs Jvector
if norm(Wnc) == 0
    Wnc = zeros(1,length(S)-1); Wnc(1) = 1;
end
if norm(Couplings) == 0
    Couplings = zeros(2,length(S)); Couplings(1,1) = 3; Couplings(:,2:end) = 1; Couplings(2,1) = 0;
end
if norm(D) == 0
    if (Nc_sim == 1)
         D = zeros(length(S),3); D(1,3) = 4;  % default anisotropy Nc
    else
         D = zeros(length(S),3); 
    end
end
if norm(J) == 0
    J = zeros(length(S),length(S));   % default J
end
if norm(Jq) == 0
    Jq = zeros(length(S),length(S));   % default Jq
end
if (Nc_sim == 1)
  if (Number_of_states_to_plot_in_spectrum == 0)
      Number_of_states_to_plot_in_spectrum = prod(2*S+1);  % here we set the default value of this variable
  end
  % construct S2 of the system without Nc:
  S2_system = sum(SxR(:,:,2:end),3)^2 + sum(SyR(:,:,2:end),3)^2 + sum(SzR(:,:,2:end),3)^2;
  % array of basic outputs:
  Er = []; spin_r = []; spin_system_r = []; sz_nc_r = []; sz_system_r = []; % array of eigenvalues and spins for each step of Jr.
  Wnc = [0, Wnc];  % weights of the couplig of Nc with the other sits. The first one is the self-coupling (Nc with Nc), that is zero.
  if (heights==0)
     Jr = linspace(Jmin,Jmax,nJ);
  else
     hr = linspace(hmin,hmax,nh);
     if (length(Jdecay)>1)
         length(hr)
         Jr = J0*exp(-Jdecay.*hr);
     else
         Jr = J0*exp(-Jdecay*hr);
     end % length(Jdecay)>1
  end
  if (norm(Jvector) > 0)
      Jr = Jvector;
  end 
  JwithNc = 0*J;
  JwithNc(1,:) = Wnc;
  % prepare spectral function for the cotunneling calculation:
  if (cotunneling == 1) % following "Theory of Single-Spin Inelastic Tunneling Spectroscopy" (J. Fernández-Rossier)
    Op = zeros(length(S2),length(S2),3);
    for isite = 1:size(SxR,3)
      Op(:,:,1) = Op(:,:,1) + Couplings(1,isite)*SxR(:,:,isite); Op(:,:,2) = Op(:,:,2) + Couplings(1,isite)*SyR(:,:,isite); Op(:,:,3) = Op(:,:,3) + Couplings(1,isite)*SzR(:,:,isite);
    end % isite
  end % if cotunneling == 1

  if (Sim_system_noNc == 1)
    Ss = S; Js = J; Jqs = Jq; Ds = D;
    Ss(1) = []; Js(1,:) = []; Js(:,1) = []; Jqs(1,:) = []; Jqs(:,1) = []; Ds(1,:) = [];
    [SxRs,SyRs,SzRs,S2s,Bs] = generate_spin_basis(Ss);
    %SxRs
    %SyRs
    %SzRs
    [Hs] = construct_spin_hamiltonian_given_coefficients(0,SxRs,SyRs,SzRs,S2s,cs2_s,Js,Jqs,Bmag,Ds);
    if (cotunneling_system_alone == 1)
      Ops = zeros(length(S2s),length(S2s),3);
      for isite = 1:size(SxRs,3)
        Ops(:,:,1) = Ops(:,:,1) + Couplings(2,isite+1)*SxRs(:,:,isite); Ops(:,:,2) = Ops(:,:,2) + Couplings(2,isite+1)*SyRs(:,:,isite); Ops(:,:,3) = Ops(:,:,3) + Couplings(2,isite+1)*SzRs(:,:,isite);
      end % isite
    end % cotunneling_system_alone
  end % Sim_system_noNc
 
  if ( Varying_Couplings > 0 )
    CouplingsR = zeros(length(Jr),length(S));
    Couplings0 = Couplings(1,:); % Initial Couplings
    Couplings1 = Couplings(3,:); % Final Couplings
    if (Varying_Couplings == 1)
      % lineal dependency
      for ii = 1:length(Jr)
        CouplingsR(ii,:) = Couplings0 + (ii - 1)*(Couplings1 - Couplings0)/(length(Jr) - 1);
      end % ii
    end % linear
    if (Varying_Couplings == 2)
      % exponential dependency
      e = 2.7172;
      a = (log(Couplings1./Couplings0))/(length(Jr) - 1); a(isnan(a)) = 0;
      for ii = 1:length(Jr)
        CouplingsR(ii,:) = Couplings0.*e.^(a*(ii-1));
      end % ii
    end % linear
    if (heights == 1)
        CouplingsR = flipud(CouplingsR);
    end
  end % if ( Varying_Couplings > 0 )

else % Nc_sim == 0
   if (cotunneling == 1)
          Op = zeros(length(S2),length(S2),3);
          for isite = 1:size(SxR,3)
           Op(:,:,1) = Op(:,:,1) + Couplings(2,isite)*SxR(:,:,isite); Op(:,:,2) = Op(:,:,2) + Couplings(2,isite)*SyR(:,:,isite); Op(:,:,3) = Op(:,:,3) + Couplings(2,isite)*SzR(:,:,isite);
           if (norm(Couplings(2,:))==0)
               Op(:,:,1) = Op(:,:,1) + SxR(:,:,isite); Op(:,:,2) = Op(:,:,2) + SyR(:,:,isite); Op(:,:,3) = Op(:,:,3) + SzR(:,:,isite);
           end % if 
          end % isite
    end % cotunneling == 1
end % if Nc_sim

end % end function set_up_variables



%%% versions of plot_spectra

function[f] = plot_spectra_h_1(Er,x,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,name_fig)
  %%% draw the branches according to their original quantum number
f=figure("visible","off"); hold on;

          for state = 1:Number_of_states_to_plot_in_spectrum
            %colormulti='r';
            %if (abs(sval(state,ij)-0.5)<=0.01); colormulti='b'; end %duplet
            %if (abs(sval(state,ij)-1.5)<=0.01); colormulti='g'; end %quadruplet
            %plot(x(ix),E(state,ij)-E(1,ij),'+','LineWidth',4,'Color',colormulti)
            %plot(x,Er(state,:),'+','LineWidth',4,'Color',colormulti)
           
            if (state<Number_of_states_to_plot_in_spectrum)
                if (  norm( Er(:,state) - Er(:,state+1) ) < 0.01   )
                else
                    if (state>1)
                        if (  norm( Er(:,state) - Er(:,state-1) ) > 0.01   )
                             name_state = strcat('Sz(Nc)=',num2str(round(sz_nc_r(1,state),3)),' S2=',num2str(round(spin_system_r(1,state),3)),' Sz=',num2str(round(sz_system_r(1,state),3)));
                        else
                            if ( sz_system_r(1,state) ~= sz_system_r(1,state-1) )
                                  name_state = strcat('Sz(Nc)=',num2str(round(sz_nc_r(1,state),3)),' S2=',num2str(round(spin_system_r(1,state),3)),' Sz=',num2str(round(sz_system_r(1,state),3)),'/',num2str(round(sz_system_r(1,state-1),3)));
                            end
                            if ( sz_nc_r(1,state) ~= sz_nc_r(1,state-1) )
                                  name_state = strcat('Sz(Nc)=',num2str(round(sz_nc_r(1,state),3)),'/',num2str(round(sz_nc_r(1,state-1),3)),' S2=',num2str(round(spin_system_r(1,state),3)),' Sz=',num2str(round(sz_system_r(1,state),3)));
                            end
                        end
                    end
                    plot(x,Er(:,state),'+','LineWidth',4,'DisplayName',name_state);
                end
            end             
           end % state

lgd = legend; lgd.FontSize = 6; lgd.Title.String = 'States'; lgd.NumColumns = 2; lgd.Location = 'northwest';

if (   abs((max(x) - min(x)) - 2) < 0.5  )
    xlabel('height (A)','FontWeight','bold','FontSize',28)
else
    xlabel('J (meV)','FontWeight','bold','FontSize',28)
end

ylabel('Eigenvalues (meV)','FontWeight','bold','FontSize',28)
set(gca, "linewidth", 4, "fontsize", 28)
saveas(f,name_fig);
end % end function plot_spectra



function[f] = plot_spectra_h_2(Er,x,Number_of_states_to_plot_in_spectrum,spin_r,spin_system_r,sz_nc_r,sz_system_r,name_fig)
  %%% draw the branches according to their original quantum number
f=figure("visible","off"); hold on;

          for state = 1:Number_of_states_to_plot_in_spectrum
           
                    c = 1*abs(sz_nc_r(:,state) - sz_nc_r(1,state));
                    d =  abs(spin_system_r(:,state) - spin_system_r(1,state));
                    %patch(x,Er(:,state),c,'EdgeColor','interp','LineWidth',4,'LineJoin','round');
                    %patch(x,Er(:,state),c,'EdgeColor','interp','LineWidth',4);
                    %patch(x,Er(:,state),c,'LineWidth',4); 
                    colorhere = [0 1 0];
                    if ( abs( spin_system_r(1,state) - 0.0 ) < 0.01 )
                        colorhere = [1 0 0];
                    end

                    if ( abs( spin_system_r(1,state) - 0.5 ) < 0.01 )
                        colorhere = [1 0 0];
                    end

                    if ( abs( spin_system_r(1,state) - 1.0 ) < 0.01 )
                        colorhere = [0.9290 0.6940 0.1250];
                    end

                    if ( abs( spin_system_r(1,state) - 1.5 ) < 0.01 )
                        colorhere = [0.9290 0.6940 0.1250];
                    end

                    if ( abs( spin_system_r(1,state) - 2.0 ) < 0.01 )
                        colorhere = [0 1 0];
                    end

                    if ( abs( spin_system_r(1,state) - 2.5 ) < 0.01 )
                        colorhere = [0 1 0];
                    end
                    

                    if (max(c)==0)
                    else
                        c = c/(max(c));
                    end

                    if (max(d)==0)
                    else
                        d = d/(max(d));
                    end
                    
                    c = 1 - c;
                    d = 1 - d;

                    for ix = 1:length(x)
                        %plot( x(ix) , Er(ix,state) )
                        %scatter(x(ix),Er(ix,state),75,"MarkerEdgeColor","b", "MarkerFaceColor",d(ix)*[0 0.7 0.7])
                        scatter(x(ix),Er(ix,state),75,"MarkerEdgeColor",d(ix)*colorhere, "MarkerFaceColor",d(ix)*colorhere)
                    end

          end % state

%lgd = legend; lgd.FontSize = 6; lgd.Title.String = 'States'; lgd.NumColumns = 2; lgd.Location = 'northwest';

if (   abs((max(x) - min(x)) - 2) < 0.5  )
    xlabel('height (A)','FontWeight','bold','FontSize',28)
else
    xlabel('J (meV)','FontWeight','bold','FontSize',28)
end

ylabel('Eigenvalues (meV)','FontWeight','bold','FontSize',28)
set(gca, "linewidth", 4, "fontsize", 28)
saveas(f,name_fig);
end % end function plot_spectra