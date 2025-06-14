module compute_lambda_and_xi
    use global_variables
    use kpoints
    implicit none
    
    contains
    
    subroutine compute_lambda_tttt(I,J,M,N,iQ_r,iQ_l,lambda1)
        
        integer,intent(in) :: I,J,M,N
        integer,intent(in) :: iQ_l,iQ_r
        integer :: c,cp,v,vp
        integer :: ik
        integer :: imQl,imQr,ikp
        double precision :: mQl(3),mQr(3),kp(3),kp2(3)
        complex(kind=16),intent(inout) :: lambda1
        complex(kind=16) :: sum1
        complex(kind=16) :: AI_ck_vkmQr,AJ_cpkp_vpkppQr,AM_ck_vpkmQl,AN_cpkp_vkppQl
        ! A_(s,v,c,k,I,Q)
        sum1 = 0
        lambda1 = 0.0
        mQl = -sys%Qpts(:,iQ_l) 
        call Qpoint_to_index(mQl,imQl)
        mQr = -sys%Qpts(:,iQ_r)
        call Qpoint_to_index(mQr,imQr)
        !mQr
        !imQr
       
        !call cpu_time(loop_start)
        do ik = 1,sys%nk !k
           kp = sys%kpts(:,ik) -sys%Qpts(:,iQ_r) - sys%Qpts(:,iQ_l)
           kp2 = kp

           call kpoint_to_index(kp,ikp)
           if(ikp == -1) then
               !print*, "Warning: kpoint index mismatch in compute_lambda_tttt"
               print*, "ik,ikp,iQ_r,iQ_l", ik,ikp,iQ_r,iQ_l
               print*, "kp", kp(1),kp(2),kp(3)
               print*, "kp2", kp2(1),kp2(2),kp2(3)
           end if
           do cp=1,sys%nc ! c'
            do vp = 1,sys%nv !v'
                do c = 1,sys%nc !c
                    do v = 1, sys%nv !v
                        
                        AI_ck_vkmQr = exciton_sys%eigenvectors_t(1,v,c,ik,I,iQ_r)
                        AJ_cpkp_vpkppQr = exciton_sys%eigenvectors_t(1,vp,cp,ikp,J,imQr)
                        AM_ck_vpkmQl = exciton_sys%eigenvectors_t(1,vp,c,ik,M,iQ_l)
                        AN_cpkp_vkppQl = exciton_sys%eigenvectors_t(1,v,cp,ikp,N,imQl)
                        sum1 = sum1 + (conjg(AN_cpkp_vkppQl) * &
                                       conjg(AM_ck_vpkmQl) * AJ_cpkp_vpkppQr * AI_ck_vkmQr)

                        
                        ! use complex conjugate for complex eigen vectors
                        ! sum2 = sum2 + (exciton_sys%eigenvectors(J,x,y)*exciton_sys%eigenvectors(I,s,t) &
                        !                *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t))
                        !print*,I,J,M,N,c,cp,v,vp,exciton_sys%eigenvectors(I,v,c),exciton_sys%eigenvectors(J,vp,cp),exciton_sys%eigenvectors(M,vp,c),exciton_sys%eigenvectors(N,v,cp)
                        
                    end do
                end do
            end do
          end do
        end do
        lambda1 = sum1
        !print*,I,J,M,N,lambda1
        
        
    end subroutine compute_lambda_tttt
   ! subroutine compute_lambda_ssss(I,J,M,N,lambda1)
   !     
   !     integer,intent(in) :: I,J,M,N
   !     integer :: c,cp,v,vp
   !     double precision,intent(inout) :: lambda1
   !     double precision :: sum1
   !     double precision :: AI_cv,AJ_cpvp,AM_cvp,AN_cpv
   !     
   !     sum1 = 0
   !     lambda1 = 0.0
   !     
   !     !call cpu_time(loop_start)
   !     do cp=1,sys%nc ! c'
   !         do vp = 1,sys%nv !v'
   !             do c = 1,sys%nc !c
   !                 do v = 1, sys%nv !v
   !                     AI_cv = exciton_sys%eigenvectors_s(I,v,c)
   !                     AJ_cpvp = exciton_sys%eigenvectors_s(J,vp,cp)
   !                     AM_cvp = exciton_sys%eigenvectors_s(M,vp,c)
   !                     AN_cpv = exciton_sys%eigenvectors_s(N,v,cp)
   !                     sum1 = sum1 + (AI_cv*AJ_cpvp *AM_cvp*AN_cpv)
   !                     
   !                     ! use complex conjugate for complex eigen vectors
   !                     ! sum2 = sum2 + (exciton_sys%eigenvectors(J,x,y)*exciton_sys%eigenvectors(I,s,t) &
   !                     !                *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t))
   !                     !print*,I,J,M,N,c,cp,v,vp,exciton_sys%eigenvectors(I,v,c),exciton_sys%eigenvectors(J,vp,cp),exciton_sys%eigenvectors(M,vp,c),exciton_sys%eigenvectors(N,v,cp)
   !                     
   !                 end do
   !             end do
   !         end do
   !     end do
   !     lambda1 = sum1
   !     !print*,I,J,M,N,lambda1
   !     
   !     
   ! end subroutine compute_lambda_ssss
!
!
   ! subroutine compute_lambda_ttss(I,J,M,N,lambda1)
   !     
   !     integer,intent(in) :: I,J,M,N
   !     integer :: c,cp,v,vp
   !     double precision,intent(inout) :: lambda1
   !     double precision :: sum1
   !     double precision :: AI_cv,AJ_cpvp,AM_cvp,AN_cpv
   !     
   !     sum1 = 0
   !     lambda1 = 0.0
   !     
   !     !call cpu_time(loop_start)
   !     do cp=1,sys%nc ! c'
   !         do vp = 1,sys%nv !v'
   !             do c = 1,sys%nc !c
   !                 do v = 1, sys%nv !v
   !                     AI_cv = exciton_sys%eigenvectors_s(I,v,c)
   !                     AJ_cpvp = exciton_sys%eigenvectors_s(J,vp,cp)
   !                     AM_cvp = exciton_sys%eigenvectors_t(M,vp,c)
   !                     AN_cpv = exciton_sys%eigenvectors_t(N,v,cp)
   !                     sum1 = sum1 + (AI_cv*AJ_cpvp *AM_cvp*AN_cpv)
   !                     
   !                     ! use complex conjugate for complex eigen vectors
   !                     ! sum2 = sum2 + (exciton_sys%eigenvectors(J,x,y)*exciton_sys%eigenvectors(I,s,t) &
   !                     !                *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t))
   !                     !print*,I,J,M,N,c,cp,v,vp,exciton_sys%eigenvectors(I,v,c),exciton_sys%eigenvectors(J,vp,cp),exciton_sys%eigenvectors(M,vp,c),exciton_sys%eigenvectors(N,v,cp)
   !                     
   !                 end do
   !             end do
   !         end do
   !     end do
   !     lambda1 = sum1
   !     !print*,I,J,M,N,lambda1
   !     
   !     
   ! end subroutine compute_lambda_ttss
!
   ! subroutine compute_lambda_sstt(I,J,M,N,lambda1)
   !     
   !     integer,intent(in) :: I,J,M,N
   !     integer :: c,cp,v,vp
   !     double precision,intent(inout) :: lambda1
   !     double precision :: sum1
   !     double precision :: AI_cv,AJ_cpvp,AM_cvp,AN_cpv
   !     
   !     sum1 = 0
   !     lambda1 = 0.0
   !     
   !     !call cpu_time(loop_start)
   !     do cp=1,sys%nc ! c'
   !         do vp = 1,sys%nv !v'
   !             do c = 1,sys%nc !c
   !                 do v = 1, sys%nv !v
   !                     AI_cv = exciton_sys%eigenvectors_t(I,v,c)
   !                     AJ_cpvp = exciton_sys%eigenvectors_t(J,vp,cp)
   !                     AM_cvp = exciton_sys%eigenvectors_s(M,vp,c)
   !                     AN_cpv = exciton_sys%eigenvectors_s(N,v,cp)
   !                     sum1 = sum1 + (AI_cv*AJ_cpvp *AM_cvp*AN_cpv)
   !                     
   !                     ! use complex conjugate for complex eigen vectors
   !                     ! sum2 = sum2 + (exciton_sys%eigenvectors(J,x,y)*exciton_sys%eigenvectors(I,s,t) &
   !                     !                *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t))
   !                     !print*,I,J,M,N,c,cp,v,vp,exciton_sys%eigenvectors(I,v,c),exciton_sys%eigenvectors(J,vp,cp),exciton_sys%eigenvectors(M,vp,c),exciton_sys%eigenvectors(N,v,cp)
   !                     
   !                 end do
   !             end do
   !         end do
   !     end do
   !     lambda1 = sum1
   !     !print*,I,J,M,N,lambda1
   !     
   !     
   ! end subroutine compute_lambda_sstt
   ! 
   ! 
   ! 
   ! subroutine compute_xi_ee(Iex,Jex,Mex,Nex,bsemat,xi_ee1_tttt,xi_ee2_tttt,xi_ee1_ssss,xi_ee2_ssss,xi_ee1_sstt,xi_ee2_sstt,xi_ee1_ttss,xi_ee2_ttss)
   !     !\xi_ee(I,J,M,N) = \sum_{i,j,c,c',v,v'}  A_{M}^{i,v}A_{N}^{j,v'}bsemat(cp,j,i,c)A_{I)^{c,v}A_{J}^{c'v'}
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,j
   !     integer :: nv
   !     double precision,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(inout) :: xi_ee1_tttt,xi_ee2_tttt,xi_ee1_ssss,xi_ee2_ssss,xi_ee1_sstt,xi_ee2_sstt,xi_ee1_ttss,xi_ee2_ttss
   !     double precision :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,AIT_cv,AIS_cv,AJT_cpvp,AJS_cpvp,AMT_iv,AMS_iv,ANT_jvp,ANS_jvp,AMT_ivp,AMS_ivp,ANT_jv,ANS_jv
   !     nv = sys%nv
   !     AIT_cv = 0.0
   !     AIS_cv = 0.0
   !     AJT_cpvp = 0.0
   !     AJS_cpvp = 0.0
   !     AMT_iv = 0.0
   !     AMS_iv = 0.0
   !     ANT_jvp = 0.0
   !     ANS_jvp = 0.0
   !     AMT_ivp = 0.0
   !     AMS_ivp = 0.0
   !     ANT_jv = 0.0
   !     ANS_jv = 0.0
   !     sum1 = 0.0
   !     sum2 = 0.0
   !     sum3 = 0.0
   !     sum4 = 0.0
   !     sum5 = 0.0
   !     sum6 = 0.0
   !     sum7 = 0.0
   !     sum8 = 0.0
   !     do i = 1,sys%nc
   !         do j = 1,sys%nc
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AIT_cv= exciton_sys%eigenvectors_t(Iex,v,c)
   !                             AIS_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJT_cpvp = exciton_sys%eigenvectors_t(Jex,vp,cp)
   !                             AJS_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AMT_iv = exciton_sys%eigenvectors_t(Mex,v,i)
   !                             AMS_iv = exciton_sys%eigenvectors_s(Mex,v,i)
   !                             ANT_jvp = exciton_sys%eigenvectors_t(Nex,vp,j)
   !                             ANS_jvp = exciton_sys%eigenvectors_s(Nex,vp,j)
   !                             AMT_ivp = exciton_sys%eigenvectors_t(Mex,vp,i)
   !                             AMS_ivp = exciton_sys%eigenvectors_s(Mex,vp,i)
   !                             ANT_jv = exciton_sys%eigenvectors_t(Nex,v,j)
   !                             ANS_jv = exciton_sys%eigenvectors_s(Nex,v,j)
   !                             sum1 = sum1 + AMT_iv*ANT_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum2 = sum2 + AMT_ivp*ANT_jv*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum3 = sum3 + AMS_iv*ANS_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum4 = sum4 + AMS_ivp*ANS_jv*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum5 = sum5 + AMS_iv*ANS_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum6 = sum6 + AMS_ivp*ANS_jv*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum7 = sum7 + AMT_iv*ANT_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum8 = sum8 + AMT_ivp*ANT_jv*bsemat(cp+nv,j+nv,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                            
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_ee1_tttt = sum1
   !     xi_ee2_tttt = sum2
   !     xi_ee1_ssss = sum3
   !     xi_ee2_ssss = sum4
   !     xi_ee1_sstt = sum5
   !     xi_ee2_sstt = sum6
   !     xi_ee1_ttss = sum7
   !     xi_ee2_ttss = sum8
   !     !xi_ee2 = sum2
   !     
   ! end subroutine compute_xi_ee
   ! 
   ! 
   ! subroutine compute_xi_hh(Iex,Jex,Mex,Nex,bsemat,xi_hh1_tttt,xi_hh2_tttt,xi_hh1_ssss,xi_hh2_ssss,xi_hh1_sstt,xi_hh2_sstt,xi_hh1_ttss,xi_hh2_ttss)
   !     !\xi_hh(I,J,M,N) = \sum_{\alpha,\beta,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{c',beta}bsemat(cp,j,i,c)A_{I)^{c,v}A_{J}^{c'v'}
   !     ! bsemat(eta,beta,alpha,gama)) = <alpha,beta|W|gama,eta>
   !     ! bsemat(beta,vp,v,alpha)) = <vvp|W|alpha,beta>
   !     
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,alpha,beta
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(inout) :: xi_hh1_tttt,xi_hh2_tttt,xi_hh1_ssss,xi_hh2_ssss,xi_hh1_sstt,xi_hh2_sstt,xi_hh1_ttss,xi_hh2_ttss
   !     double precision :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,AIT_cv,AIS_cv,AJT_cpvp,AJS_cpvp,AMT_calpha,AMS_calpha,ANT_cpbeta,ANS_cpbeta,AMT_cbeta,AMS_cbeta,ANT_cpalpha,ANS_cpalpha
   !     sum1 = 0
   !     sum2 = 0
   !     sum3 = 0
   !     sum4 = 0
   !     sum5 = 0
   !     sum6 = 0
   !     sum7 = 0
   !     sum8 = 0
   !     xi_hh1_tttt=0.0
   !     xi_hh2_tttt=0.0
   !     xi_hh1_ssss=0.0
   !     xi_hh2_ssss=0.0
   !     xi_hh1_sstt=0.0
   !     xi_hh2_sstt=0.0
   !     xi_hh1_ttss=0.0
   !     xi_hh2_ttss=0.0
   !     do alpha = 1,sys%nv
   !         do beta = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AIT_cv= exciton_sys%eigenvectors_t(Iex,v,c)
   !                             AIS_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJT_cpvp = exciton_sys%eigenvectors_t(Jex,vp,cp)
   !                             AJS_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AMT_calpha = exciton_sys%eigenvectors_t(Mex,alpha,c)
   !                             AMS_calpha = exciton_sys%eigenvectors_s(Mex,alpha,c)
   !                             ANT_cpbeta = exciton_sys%eigenvectors_t(Nex,beta,cp)
   !                             ANS_cpbeta = exciton_sys%eigenvectors_s(Nex,beta,cp)
   !                             AMT_cbeta = exciton_sys%eigenvectors_t(Mex,beta,c)
   !                             AMS_cbeta = exciton_sys%eigenvectors_s(Mex,beta,c)
   !                             ANT_cpalpha = exciton_sys%eigenvectors_t(Nex,alpha,cp)
   !                             ANS_cpalpha = exciton_sys%eigenvectors_s(Nex,alpha,cp) 
   !                             sum1 = sum1 + AMT_calpha*ANT_cpbeta*bsemat(beta,vp,v,alpha)*AIT_cv*AJT_cpvp
   !                             sum2 = sum2 + AMT_cbeta*ANT_cpalpha*bsemat(beta,vp,v,alpha)*AIT_cv*AJT_cpvp
   !                             sum3 = sum3 + AMS_calpha*ANS_cpbeta*bsemat(beta,vp,v,alpha)*AIS_cv*AJS_cpvp
   !                             sum4 = sum4 + AMS_cbeta*ANS_cpalpha*bsemat(beta,vp,v,alpha)*AIS_cv*AJS_cpvp
   !                             sum5 = sum5 + AMS_calpha*ANS_cpbeta*bsemat(beta,vp,v,alpha)*AIT_cv*AJT_cpvp
   !                             sum6 = sum6 + AMS_cbeta*ANS_cpalpha*bsemat(beta,vp,v,alpha)*AIT_cv*AJT_cpvp
   !                             sum7 = sum7 + AMT_calpha*ANT_cpbeta*bsemat(beta,vp,v,alpha)*AIS_cv*AJS_cpvp
   !                             sum8 = sum8 + AMT_cbeta*ANT_cpalpha*bsemat(beta,vp,v,alpha)*AIS_cv*AJS_cpvp                            
   !                             
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_hh1_tttt = sum1
   !     xi_hh2_tttt = sum2
   !     xi_hh1_ssss = sum3
   !     xi_hh2_ssss = sum4
   !     xi_hh1_sstt = sum5
   !     xi_hh2_sstt = sum6
   !     xi_hh1_ttss = sum7
   !     xi_hh2_ttss = sum8
   !     
   ! end subroutine compute_xi_hh
   ! 
   ! subroutine compute_xi_in_ehd(Iex,Jex,Mex,Nex,bsemat,xi_in_ehd1_tttt,xi_in_ehd2_tttt,xi_in_ehd1_ssss,xi_in_ehd2_ssss,xi_in_ehd1_sstt,xi_in_ehd2_sstt,xi_in_ehd1_ttss,xi_in_ehd2_ttss)
   !     !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !\xi_out_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{i,v'}bsemat(alpha,v,i+nv,c'+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !bsemat(alpha,v',i+nv,c+nv) = <iv'|W|c,alpha>
   !     !bsemat(alpha,v,i+nv,c'+nv)= <iv|W|c',alpha>
   !     
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(out) :: xi_in_ehd1_tttt,xi_in_ehd2_tttt,xi_in_ehd1_ssss,xi_in_ehd2_ssss,xi_in_ehd1_sstt,xi_in_ehd2_sstt,xi_in_ehd1_ttss,xi_in_ehd2_ttss
   !     double precision :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,AIT_cv,AIS_cv,AJT_cpvp,AJS_cpvp,AMT_iv,AMS_iv,ANT_cpalpha,ANS_cpalpha,AMT_ialpha,AMS_ialpha,ANT_cpv,ANS_cpv
   !     sum1 = 0
   !     sum2 = 0
   !     sum3 = 0
   !     sum4 = 0
   !     sum5 = 0
   !     sum6 = 0
   !     sum7 = 0
   !     sum8 = 0
   !     xi_in_ehd1_tttt= 0.0
   !     xi_in_ehd2_tttt= 0.0
   !     xi_in_ehd1_ssss= 0.0
   !     xi_in_ehd2_ssss= 0.0
   !     xi_in_ehd1_sstt= 0.0
   !     xi_in_ehd2_sstt= 0.0
   !     xi_in_ehd1_ttss= 0.0
   !     xi_in_ehd2_ttss= 0.0
   !  
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AIT_cv= exciton_sys%eigenvectors_t(Iex,v,c)
   !                             AIS_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJT_cpvp = exciton_sys%eigenvectors_t(Jex,vp,cp)
   !                             AJS_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AMT_iv = exciton_sys%eigenvectors_t(Mex,v,i)
   !                             AMS_iv = exciton_sys%eigenvectors_s(Mex,v,i)
   !                             ANT_cpalpha = exciton_sys%eigenvectors_t(Nex,alpha,cp)
   !                             ANS_cpalpha = exciton_sys%eigenvectors_s(Nex,alpha,cp)
   !                             AMT_ialpha = exciton_sys%eigenvectors_t(Mex,alpha,i)
   !                             AMS_ialpha = exciton_sys%eigenvectors_s(Mex,alpha,i)
   !                             ANT_cpv = exciton_sys%eigenvectors_t(Nex,v,cp)
   !                             ANS_cpv = exciton_sys%eigenvectors_s(Nex,v,cp)
   !                             sum1 = sum1 + AMT_iv*ANT_cpalpha*bsemat(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum2 = sum2 + AMT_ialpha*ANT_cpv*bsemat(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum3 = sum3 + AMS_iv*ANS_cpalpha*bsemat(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum4 = sum4 + AMS_ialpha*ANS_cpv*bsemat(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum5 = sum5 + AMS_iv*ANS_cpalpha*bsemat(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum6 = sum6 + AMS_ialpha*ANS_cpv*bsemat(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum7 = sum7 + AMT_iv*ANT_cpalpha*bsemat(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum8 = sum8 + AMT_ialpha*ANT_cpv*bsemat(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             !sum2 = sum2 + AM_iv*AN_calpha*bsemat(alpha,vp,i+nv,c+nv)*AJ_cv*AI_cpvp
   !                             !sum3 = sum3 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
   !                             !sum4 = sum4 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AJ_cv*AI_cpvp
   !                             !print*,sum1,"sum1"
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_in_ehd1_tttt = sum1
   !     xi_in_ehd2_tttt = sum2
   !     xi_in_ehd1_ssss = sum3
   !     xi_in_ehd2_ssss = sum4
   !     xi_in_ehd1_sstt = sum5
   !     xi_in_ehd2_sstt = sum6
   !     xi_in_ehd1_ttss = sum7
   !     xi_in_ehd2_ttss = sum8
   !    ! print*,sum1
   !     !xi_in_ehd2 = sum2
   !     !xi_out_ehd1 = sum3
   !     !xi_out_ehd2 = sum4
   ! end subroutine compute_xi_in_ehd
   ! 
   ! subroutine compute_xi_out_ehd(Iex,Jex,Mex,Nex,bsemat,xi_out_ehd1_tttt,xi_out_ehd2_tttt,xi_out_ehd1_ssss,xi_out_ehd2_ssss,xi_out_ehd1_sstt,xi_out_ehd2_sstt,xi_out_ehd1_ttss,xi_out_ehd2_ttss)
   !     !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !\xi_out_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{i,v'}bsemat(alpha,v,i+nv,c'+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !bsemat(alpha,v',i+nv,c+nv) = <iv'|W|c,alpha>
   !     !bsemat(alpha,v,i+nv,c'+nv)= <iv|W|c',alpha>
   !     
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(out) :: xi_out_ehd1_tttt,xi_out_ehd2_tttt,xi_out_ehd1_ssss,xi_out_ehd2_ssss,xi_out_ehd1_sstt,xi_out_ehd2_sstt,xi_out_ehd1_ttss,xi_out_ehd2_ttss !,xi_in_ehd2,xi_out_ehd1,xi_out_ehd2
   !     double precision :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,AIT_cv,AIS_cv,AJT_cpvp,AJS_cpvp,AMT_calpha,AMS_calpha,ANT_ivp,ANS_ivp,AMT_cvp,AMS_cvp,ANT_ialpha,ANS_ialpha
   !     sum1 = 0
   !     sum2 = 0
   !     sum3 = 0
   !     sum4 = 0
   !     sum5 = 0
   !     sum6 = 0
   !     sum7 = 0
   !     sum8 = 0
   !     xi_out_ehd1_tttt = 0.0
   !     xi_out_ehd2_tttt = 0.0
   !     xi_out_ehd1_ssss = 0.0  
   !     xi_out_ehd2_ssss = 0.0
   !     xi_out_ehd1_sstt = 0.0
   !     xi_out_ehd2_sstt = 0.0  
   !     xi_out_ehd1_ttss = 0.0
   !     xi_out_ehd2_ttss = 0.0  
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AIT_cv= exciton_sys%eigenvectors_t(Iex,v,c)
   !                             AIS_cv= exciton_sys%eigenvectors_s(Iex,v,c) 
   !                             AJT_cpvp = exciton_sys%eigenvectors_t(Jex,vp,cp)
   !                             AJS_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AMT_calpha = exciton_sys%eigenvectors_t(Mex,alpha,c)
   !                             AMS_calpha = exciton_sys%eigenvectors_s(Mex,alpha,c)
   !                             ANT_ivp = exciton_sys%eigenvectors_t(Nex,vp,i)
   !                             ANS_ivp = exciton_sys%eigenvectors_s(Nex,vp,i)
   !                             AMT_cvp = exciton_sys%eigenvectors_t(Mex,vp,c)
   !                             AMS_cvp = exciton_sys%eigenvectors_s(Mex,vp,c)
   !                             ANT_ialpha = exciton_sys%eigenvectors_t(Nex,alpha,i)
   !                             ANS_ialpha = exciton_sys%eigenvectors_s(Nex,alpha,i)
   !                             sum1 = sum1 + AMT_calpha*ANT_ivp*bsemat(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum2 = sum2 + AMT_cvp*ANT_ialpha*bsemat(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum3 = sum3 + AMS_calpha*ANS_ivp*bsemat(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             sum4 = sum4 + AMS_cvp*ANS_ialpha*bsemat(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             sum5 = sum5 + AMS_calpha*ANS_ivp*bsemat(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum6 = sum6 + AMS_cvp*ANS_ialpha*bsemat(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum7 = sum7 + AMT_calpha*ANT_ivp*bsemat(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             sum8 = sum8 + AMT_cvp*ANT_ialpha*bsemat(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             !sum4 = sum4 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AJ_cv*AI_cpvp
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_out_ehd1_tttt = sum1
   !     xi_out_ehd2_tttt = sum2 
   !     xi_out_ehd1_ssss = sum3 
   !     xi_out_ehd2_ssss = sum4
   !     xi_out_ehd1_sstt = sum5 
   !     xi_out_ehd2_sstt = sum6
   !     xi_out_ehd1_ttss = sum7 
   !     xi_out_ehd2_ttss = sum8
   !     !xi_out_ehd2 = sum4
   ! end subroutine compute_xi_out_ehd
!
   ! 
!
   ! !singlet xi direct elements
!
   ! subroutine compute_xi_ee_ssss(Iex,Jex,Mex,Nex,bsemat,xi_ee1,xi_ee2)
   !     !\xi_ee(I,J,M,N) = \sum_{i,j,c,c',v,v'}  A_{M}^{i,v}A_{N}^{j,v'}bsemat(cp,j,i,c)A_{I)^{c,v}A_{J}^{c'v'}
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,j
   !     integer :: nv
   !     double precision,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(inout) :: xi_ee1,xi_ee2
   !     double precision :: sum1,sum2,AI_cv,AJ_cpvp,AM_iv,AN_jvp,AM_ivp,AN_jv
   !     nv = sys%nv
   !     AI_cv = 0.0
   !     AJ_cpvp = 0.0
   !     AM_iv = 0.0
   !     AN_jvp = 0.0
   !     AM_ivp = 0.0
   !     AN_jv = 0.0
   !     sum1 = 0.0
   !     sum2 = 0.0
   !     do i = 1,sys%nc
   !         do j = 1,sys%nc
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AI_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJ_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AM_iv = exciton_sys%eigenvectors_s(Mex,v,i)
   !                             AN_jvp = exciton_sys%eigenvectors_s(Nex,vp,j)
   !                             AM_ivp = exciton_sys%eigenvectors_s(Mex,vp,i)
   !                             AN_jv = exciton_sys%eigenvectors_s(Nex,v,j)
   !                             sum1 = sum1 + AM_iv*AN_jvp*bsemat(cp+nv,j+nv,i+nv,c+nv)*AI_cv*AJ_cpvp
   !                             sum2 = sum2 + AM_ivp*AN_jv*bsemat(cp+nv,j+nv,i+nv,c+nv)*AI_cv*AJ_cpvp
   !                            
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_ee1 = sum1
   !     xi_ee2 = sum2
   !     !xi_ee2 = sum2
   !     
   ! end subroutine compute_xi_ee_ssss
   ! 
   ! 
   ! subroutine compute_xi_hh_ssss(Iex,Jex,Mex,Nex,bsemat,xi_hh1,xi_hh2)
   !     !\xi_hh(I,J,M,N) = \sum_{\alpha,\beta,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{c',beta}bsemat(cp,j,i,c)A_{I)^{c,v}A_{J}^{c'v'}
   !     ! bsemat(eta,beta,alpha,gama)) = <alpha,beta|W|gama,eta>
   !     ! bsemat(beta,vp,v,alpha)) = <vvp|W|alpha,beta>
   !     
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,alpha,beta
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(inout) :: xi_hh1,xi_hh2
   !     double precision :: sum1,sum2,AI_cv,AJ_cpvp,AM_calpha,AN_cpbeta,AJ_cv,AI_cpvp,AM_cbeta,AN_cpalpha
   !     sum1 = 0
   !     sum2 = 0
   !     do alpha = 1,sys%nv
   !         do beta = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AI_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJ_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AM_calpha = exciton_sys%eigenvectors_s(Mex,alpha,c)
   !                             AN_cpbeta = exciton_sys%eigenvectors_s(Nex,beta,cp)
   !                             AM_cbeta = exciton_sys%eigenvectors_s(Mex,beta,c)
   !                             AN_cpalpha = exciton_sys%eigenvectors_s(Nex,alpha,cp) 
   !                             sum1 = sum1 + AM_calpha*AN_cpbeta*bsemat(beta,vp,v,alpha)*AI_cv*AJ_cpvp
   !                             sum2 = sum2 + AM_cbeta*AN_cpalpha*bsemat(beta,vp,v,alpha)*AI_cv*AJ_cpvp
   !                             
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_hh1 = sum1
   !     xi_hh2 = sum2
   !     
   ! end subroutine compute_xi_hh_ssss
   ! 
   ! subroutine compute_xi_in_ehd_ssss(Iex,Jex,Mex,Nex,bsemat,xi_in_ehd1,xi_in_ehd2)
   !     !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{cp,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !\xi_out_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{cp,alpha}A_{N}^{i,v'}bsemat(alpha,v,i+nv,c'+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !bsemat(alpha,v',i+nv,c+nv) = <iv'|W|c,alpha>
   !     !bsemat(alpha,v,i+nv,c'+nv)= <iv|W|c',alpha>
   !     
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(out) :: xi_in_ehd1,xi_in_ehd2 !,xi_in_ehd2,xi_out_ehd1,xi_out_ehd2
   !     double precision :: sum1,sum2,sum3,sum4,AI_cv,AJ_cpvp,AM_iv,AN_cpalpha,AJ_cv,AI_cpvp,AM_ialpha,AN_cpv
   !     sum1 = 0
   !     sum2 = 0
   !     xi_in_ehd1 = 0.0
   !     xi_in_ehd2 = 0.0
   !  
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AI_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJ_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AM_iv = exciton_sys%eigenvectors_s(Mex,v,i)
   !                             AN_cpalpha = exciton_sys%eigenvectors_s(Nex,alpha,cp)
   !                             AM_ialpha = exciton_sys%eigenvectors_s(Mex,alpha,i)
   !                             AN_cpv = exciton_sys%eigenvectors_s(Nex,v,cp)
   !                             sum1 = sum1 + AM_iv*AN_cpalpha*bsemat(alpha,vp,i+nv,c+nv)*AI_cv*AJ_cpvp
   !                             sum2 = sum2 + AM_ialpha*AN_cpv*bsemat(alpha,vp,i+nv,c+nv)*AI_cv*AJ_cpvp
   !                             !sum2 = sum2 + AM_iv*AN_calpha*bsemat(alpha,vp,i+nv,c+nv)*AJ_cv*AI_cpvp
   !                             !sum3 = sum3 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
   !                             !sum4 = sum4 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AJ_cv*AI_cpvp
   !                             !print*,sum1,"sum1"
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_in_ehd1 = sum1
   !     xi_in_ehd2 = sum2
   !    ! print*,sum1
   !     !xi_in_ehd2 = sum2
   !     !xi_out_ehd1 = sum3
   !     !xi_out_ehd2 = sum4
   ! end subroutine compute_xi_in_ehd_ssss
   ! 
   ! subroutine compute_xi_out_ehd_ssss(Iex,Jex,Mex,Nex,bsemat,xi_out_ehd1,xi_out_ehd2)
   !     !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !\xi_out_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{c,alpha}A_{N}^{i,v'}bsemat(alpha,v,i+nv,c'+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     !bsemat(alpha,v',i+nv,c+nv) = <iv'|W|c,alpha>
   !     !bsemat(alpha,v,i+nv,c'+nv)= <iv|W|c',alpha>
   !     
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat
   !     double precision,intent(out) :: xi_out_ehd1,xi_out_ehd2 !,xi_in_ehd2,xi_out_ehd1,xi_out_ehd2
   !     double precision :: sum1,sum2,sum3,sum4,AI_cv,AJ_cpvp,AM_calpha,AN_ivp,AM_cvp,AN_ialpha
   !     sum1 = 0
   !     sum2 = 0
   !     xi_out_ehd1 = 0.0
   !     xi_out_ehd2 = 0.0
   !   
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             !bsemat(n+nv,j+nv,i+nv,m+nv) = <ij|W|mn>
   !                             !bsemat()   = <ij|W|ccp>
   !                             AI_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJ_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AM_calpha = exciton_sys%eigenvectors_s(Mex,alpha,c)
   !                             AN_ivp = exciton_sys%eigenvectors_s(Nex,vp,i)
   !                             AM_cvp = exciton_sys%eigenvectors_s(Mex,vp,c)
   !                             AN_ialpha = exciton_sys%eigenvectors_s(Nex,alpha,i)
   !                             sum1 = sum1 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
   !                             sum2 = sum2 + AM_cvp*AN_ialpha*bsemat(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
   !                             !sum4 = sum4 + AM_calpha*AN_ivp*bsemat(alpha,v,i+nv,cp+nv)*AJ_cv*AI_cpvp
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_out_ehd1 = sum1
   !     xi_out_ehd2 = sum2 
   !     !xi_out_ehd2 = sum4
   ! end subroutine compute_xi_out_ehd_ssss
   ! 
   ! subroutine compute_xi_in_ehx(I_ex,J_ex,M_ex,N_ex,bsemat_x,xi_in_ehx1_tttt,xi_in_ehx2_tttt,xi_in_ehx1_ssss,xi_in_ehx2_ssss,xi_in_ehx1_sstt,xi_in_ehx2_sstt,xi_in_ehx1_ttss,xi_in_ehx2_ttss)
   !     !xi_in_ehx(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{cp,alpha}bsemat_x(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !      !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{cp,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     integer,intent(in) :: I_ex,J_ex,M_ex,N_ex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat_x
   !     double precision,intent(out) :: xi_in_ehx1_tttt,xi_in_ehx2_tttt,xi_in_ehx1_ssss,xi_in_ehx2_ssss,xi_in_ehx1_sstt,xi_in_ehx2_sstt,xi_in_ehx1_ttss,xi_in_ehx2_ttss
   !     double precision :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,AIT_cv,AIS_cv,AJT_cpvp,AJS_cpvp,AMT_iv,AMS_iv,ANT_cpalpha,ANS_cpalpha,AMT_ialpha,AMS_ialpha,ANT_cpv,ANS_cpv
   !     sum1 = 0
   !     sum2 = 0
   !     sum3 = 0
   !     sum4 = 0
   !     sum5 = 0
   !     sum6 = 0
   !     sum7 = 0
   !     sum8 = 0
   !     xi_in_ehx1_tttt = 0.0
   !     xi_in_ehx2_tttt = 0.0
   !     xi_in_ehx1_ssss = 0.0
   !     xi_in_ehx2_ssss = 0.0
   !     xi_in_ehx1_sstt = 0.0
   !     xi_in_ehx2_sstt = 0.0
   !     xi_in_ehx1_ttss = 0.0
   !     xi_in_ehx2_ttss = 0.0
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             AIT_cv= exciton_sys%eigenvectors_t(I_ex,v,c)
   !                             AIS_cv= exciton_sys%eigenvectors_s(I_ex,v,c)
   !                             AJT_cpvp = exciton_sys%eigenvectors_t(J_ex,vp,cp)
   !                             AJS_cpvp = exciton_sys%eigenvectors_s(J_ex,vp,cp)
   !                             AMT_iv = exciton_sys%eigenvectors_t(M_ex,v,i)
   !                             AMS_iv = exciton_sys%eigenvectors_s(M_ex,v,i)
   !                             ANT_cpalpha = exciton_sys%eigenvectors_t(N_ex,alpha,cp)
   !                             ANS_cpalpha = exciton_sys%eigenvectors_s(N_ex,alpha,cp)
   !                             AMT_ialpha = exciton_sys%eigenvectors_t(M_ex,alpha,i)
   !                             AMS_ialpha = exciton_sys%eigenvectors_s(M_ex,alpha,i)
   !                             ANT_cpv = exciton_sys%eigenvectors_t(N_ex,v,cp)
   !                             ANS_cpv = exciton_sys%eigenvectors_s(N_ex,v,cp)
   !                             sum1 = sum1 + AMT_iv*ANT_cpalpha*bsemat_x(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum2 = sum2 + AMT_ialpha*ANT_cpv*bsemat_x(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum3 = sum3 + AMS_iv*ANS_cpalpha*bsemat_x(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum4 = sum4 + AMS_ialpha*ANS_cpv*bsemat_x(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum5 = sum5 + AMS_iv*ANS_cpalpha*bsemat_x(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum6 = sum6 + AMS_ialpha*ANS_cpv*bsemat_x(alpha,vp,i+nv,c+nv)*AIT_cv*AJT_cpvp
   !                             sum7 = sum7 + AMT_iv*ANT_cpalpha*bsemat_x(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             sum8 = sum8 + AMT_ialpha*ANT_cpv*bsemat_x(alpha,vp,i+nv,c+nv)*AIS_cv*AJS_cpvp
   !                             if(I_ex==4 .and. J_ex==1 .and. M_ex==4 .and. N_ex==1) then
   !                                 !print*,"AAAAb",AM_iv,AN_cpalpha,bsemat_x(alpha,vp,i+nv,c+nv)*AI_cv,alpha,vp,i,c
   !                             end if
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do 
   !     end do
   !    xi_in_ehx1_tttt = sum1
   !    xi_in_ehx2_tttt = sum2
   !    xi_in_ehx1_ssss = sum3
   !    xi_in_ehx2_ssss = sum4
   !    xi_in_ehx1_sstt = sum5
   !    xi_in_ehx2_sstt = sum6
   !    xi_in_ehx1_ttss = sum7
   !    xi_in_ehx2_ttss = sum8
!
   ! end subroutine compute_xi_in_ehx
!
!
   ! subroutine compute_xi_out_ehx(Iex,Jex,Mex,Nex,bsemat_x,xi_out_ehx1_tttt,xi_out_ehx2_tttt,xi_out_ehx1_ssss,xi_out_ehx2_ssss,xi_out_ehx1_sstt,xi_out_ehx2_sstt,xi_out_ehx1_ttss,xi_out_ehx2_ttss)
   !     !xi_in_ehx(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat_x(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
!
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat_x
   !     double precision,intent(inout) :: xi_out_ehx1_tttt,xi_out_ehx2_tttt,xi_out_ehx1_ssss,xi_out_ehx2_ssss,xi_out_ehx1_sstt,xi_out_ehx2_sstt,xi_out_ehx1_ttss,xi_out_ehx2_ttss
   !     double precision :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,AIT_cv,AIS_cv,AJT_cpvp,AJS_cpvp,AMT_calpha,AMS_calpha,ANT_ivp,ANS_ivp,AMT_cvp,AMS_cvp,ANT_ialpha,ANS_ialpha
   !     
   !     sum1 = 0
   !     sum2 = 0
   !     sum3 = 0
   !     sum4 = 0
   !     sum5 = 0
   !     sum6 = 0
   !     sum7 = 0
   !     sum8 = 0
   !     xi_out_ehx1_tttt = 0.0
   !     xi_out_ehx2_tttt = 0.0
   !     xi_out_ehx1_ssss = 0.0
   !     xi_out_ehx2_ssss = 0.0
   !     xi_out_ehx1_sstt = 0.0
   !     xi_out_ehx2_sstt = 0.0
   !     xi_out_ehx1_ttss = 0.0
   !     xi_out_ehx2_ttss = 0.0
!
!
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             AIT_cv= exciton_sys%eigenvectors_t(Iex,v,c)
   !                             AIS_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJT_cpvp = exciton_sys%eigenvectors_t(Jex,vp,cp)
   !                             AJS_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AMT_calpha = exciton_sys%eigenvectors_t(Mex,alpha,c)
   !                             AMS_calpha = exciton_sys%eigenvectors_s(Mex,alpha,c)
   !                             ANT_ivp = exciton_sys%eigenvectors_t(Nex,vp,i)
   !                             ANS_ivp = exciton_sys%eigenvectors_s(Nex,vp,i)
   !                             AMT_cvp = exciton_sys%eigenvectors_t(Mex,vp,c)
   !                             AMS_cvp = exciton_sys%eigenvectors_s(Mex,vp,c)
   !                             ANT_ialpha = exciton_sys%eigenvectors_t(Nex,alpha,i)
   !                             ANS_ialpha = exciton_sys%eigenvectors_s(Nex,alpha,i)
   !                             sum1 = sum1 + AMT_calpha*ANT_ivp*bsemat_x(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum2 = sum2 + AMT_cvp*ANT_ialpha*bsemat_x(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum3 = sum3 + AMS_calpha*ANS_ivp*bsemat_x(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             sum4 = sum4 + AMS_cvp*ANS_ialpha*bsemat_x(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             sum5 = sum5 + AMS_calpha*ANS_ivp*bsemat_x(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum6 = sum6 + AMS_cvp*ANS_ialpha*bsemat_x(alpha,v,i+nv,cp+nv)*AIT_cv*AJT_cpvp
   !                             sum7 = sum7 + AMT_calpha*ANT_ivp*bsemat_x(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                             sum8 = sum8 + AMT_cvp*ANT_ialpha*bsemat_x(alpha,v,i+nv,cp+nv)*AIS_cv*AJS_cpvp
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_out_ehx1_tttt = sum1
   !     xi_out_ehx2_tttt = sum2
   !     xi_out_ehx1_ssss = sum3
   !     xi_out_ehx2_ssss = sum4
   !     xi_out_ehx1_sstt = sum5
   !     xi_out_ehx2_sstt = sum6
   !     xi_out_ehx1_ttss = sum7
   !     xi_out_ehx2_ttss = sum8
!
   ! end subroutine compute_xi_out_ehx
!
   ! subroutine compute_xi_in_ehx_ssss(I_ex,J_ex,M_ex,N_ex,bsemat_x,xi_in_ehx1,xi_in_ehx2)
   !     !xi_in_ehx(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{cp,alpha}bsemat_x(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !      !\xi_in_ehd(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{cp,alpha}bsemat(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
   !     integer,intent(in) :: I_ex,J_ex,M_ex,N_ex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat_x
   !     double precision,intent(out) :: xi_in_ehx1,xi_in_ehx2
   !     double precision :: sum1,sum2,AI_cv,AJ_cpvp,AM_iv,AN_cpalpha,AM_ialpha,AN_cpv
   !     sum1 = 0
   !     sum2 = 0
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             AI_cv= exciton_sys%eigenvectors_s(I_ex,v,c)
   !                             AJ_cpvp = exciton_sys%eigenvectors_s(J_ex,vp,cp)
   !                             AM_iv = exciton_sys%eigenvectors_s(M_ex,v,i)
   !                             AN_cpalpha = exciton_sys%eigenvectors_s(N_ex,alpha,cp)
   !                             AM_ialpha = exciton_sys%eigenvectors_s(M_ex,alpha,i)
   !                             AN_cpv = exciton_sys%eigenvectors_s(N_ex,v,cp)
   !                             sum1 = sum1 + AM_iv*AN_cpalpha*bsemat_x(alpha,vp,i+nv,c+nv)*AI_cv*AJ_cpvp
   !                             sum2 = sum2 + AM_ialpha*AN_cpv*bsemat_x(alpha,vp,i+nv,c+nv)*AI_cv*AJ_cpvp
   !                             if(I_ex==4 .and. J_ex==1 .and. M_ex==4 .and. N_ex==1) then
   !                                 !print*,"AAAAb",AM_iv,AN_cpalpha,bsemat_x(alpha,vp,i+nv,c+nv)*AI_cv,alpha,vp,i,c
   !                             end if
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do 
   !     end do
!
   !     xi_in_ehx1 = sum1
   !     xi_in_ehx2 = sum2
!
   ! end subroutine compute_xi_in_ehx_ssss
!
!
   ! subroutine compute_xi_out_ehx_ssss(Iex,Jex,Mex,Nex,bsemat_x,xi_out_ehx1,xi_out_ehx2)
   !     !xi_in_ehx(I,J,M,N) = \sum_{i,\alpha,c,c',v,v'}  A_{M}^{i,v}A_{N}^{c,alpha}bsemat_x(alpha,v',i+nv,c+nv)A_{I)^{c,v}A_{J}^{c'v'}
!
   !     integer,intent(in) :: Iex,Jex,Mex,Nex
   !     integer :: c,v,cp,vp,i,alpha,nv
   !     double precision ,intent(in),dimension(sys%nb,sys%nb,sys%nb,sys%nb) :: bsemat_x
   !     double precision,intent(inout) :: xi_out_ehx1,xi_out_ehx2
   !     double precision :: sum1,sum2,sum3,sum4,AI_cv,AJ_cpvp,AM_calpha,AN_ivp,AM_cvp,AN_ialpha
   !     
   !     sum1 = 0
   !     sum2 = 0
   !     nv = sys%nv
   !     do i = 1,sys%nc  !
   !         do alpha = 1,sys%nv
   !             do c = 1,sys%nc
   !                 do cp = 1,sys%nc
   !                     do v = 1,sys%nv
   !                         do vp = 1,sys%nv
   !                             AI_cv= exciton_sys%eigenvectors_s(Iex,v,c)
   !                             AJ_cpvp = exciton_sys%eigenvectors_s(Jex,vp,cp)
   !                             AM_calpha = exciton_sys%eigenvectors_s(Mex,alpha,c)
   !                             AN_ivp = exciton_sys%eigenvectors_s(Nex,vp,i)
   !                             AM_cvp = exciton_sys%eigenvectors_s(Mex,vp,c)
   !                             AN_ialpha = exciton_sys%eigenvectors_s(Nex,alpha,i)
   !                             sum1 = sum1 + AM_calpha*AN_ivp*bsemat_x(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
   !                             sum2 = sum2 + AM_cvp*AN_ialpha*bsemat_x(alpha,v,i+nv,cp+nv)*AI_cv*AJ_cpvp
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   !     xi_out_ehx1 = sum1
   !     xi_out_ehx2 = sum2
!
   ! end subroutine compute_xi_out_ehx_ssss
!
!
   ! 
!
   !      
   ! 
   ! 
   !! subroutine compute_xi_in()
   !!     
   !!     implicit none
   !!     integer :: I,J,M,N,R,S
   !!     
   !!     do M = 1,sys%nex
   !!         do N = 1,sys%nex
   !!             do I= 1,sys%nex
   !!                 do J = 1,sys%nex
   !!                     do R = 1,sys%nex
   !!                         do S = 1,sys%nex
   !!                             ham%xi_in(I,J,M,N) = ham%xi_in(I,J,M,N) + ham%lambda(R,S,M,N)*ham%xi(I,J,R,S)
   !!                             !
   !!                         end do
   !!                     end do
   !!                    !print*,"xi_in",I,J,M,N,ham%xi_in(I,J,M,N)
   !!                 end do
   !!             end do
   !!         end do
   !!     end do
   !!     
   !!     
   !!     
   ! end subroutine compute_xi_in
   ! 
   ! subroutine compute_xi_in_min_xi_dir_x_lam_eph()
   !     
   !     implicit none
   !     integer :: I,J,M,N,R,S
   !     print*,"hello2"
   !     
   !     do I = 1,sys%nex
   !         do J = 1,sys%nex
   !             do M = 1,sys%nex
   !                 do N = 1,sys%nex
   !                     do R = 1,sys%nex
   !                         do S = 1,sys%nex
   !                             
   !                             ham%xi_in_min_xi_dir_x_lam_eph(I,J,M,N) = ham%xi_in_min_xi_dir_x_lam_eph(I,J,M,N) &
   !                             + ((ham%xi_in(R,S,M,N) - ham%xi(R,S,M,N))*(ham%lambda_e(I,J,R,S) + ham%lambda(I,J,R,S)))
   !                         end do
   !                     end do
   !                     !print*,I,J,M,N,ham%xi_in_min_xi_dir_x_lam_eph(I,J,M,N)
   !                 end do
   !             end do
   !         end do
   !     end do
   !     
   ! end subroutine compute_xi_in_min_xi_dir_x_lam_eph
   ! 
   ! subroutine compute_lambda_h_ERS_lambda_eph()
   !     
   !     implicit none
   !     integer :: I,J,M,N,R,S,c2
   !     
   !     do I = 1,sys%nex
   !         do J = 1,sys%nex
   !             do M = 1,sys%nex
   !                 do N = 1,sys%nex
   !                     c2 = 0
   !                     do R = 1,sys%nex
   !                         do S = 1,sys%nex
   !                             ham%lambda_h_ERS_lambda_eph(I,J,M,N) = ham%lambda_h_ERS_lambda_eph(I,J,M,N) &
   !                                   + (ham%lambda(R,S,M,N)*(exciton_sys%eigenvalues(R)+exciton_sys%eigenvalues(S)) &
   !                                   * (ham%lambda_e(I,J,R,S)+ham%lambda(I,J,R,S)))   
   !                             c2 = c2 + 1
   !                             !print*,"first gather",c2,(ham%lambda(R,S,M,N)*(exciton_sys%eigenvalues(R)+exciton_sys%eigenvalues(S)))                             
   !                         end do
   !                     end do
   !                     !print*,"lambda_h_ERS_lambda_eph",I,J,M,N,ham%lambda_h_ERS_lambda_eph(I,J,M,N)
   !                 end do
   !             end do
   !         end do
   !     end do
!
   ! end subroutine compute_lambda_h_ERS_lambda_eph
   !     
   !     !subroutine compute_xi_out()
   !     !   implicit none
   !     !   integer :: I,J,M,N,R,S
   !     !commented out the following block of code as the subroutine does not seem to be used anywhere
   !     !subroutine compute_xi_out()
   !     !    implicit none
   !     !    integer :: I,J,M,N,R,S
   !     
   !     !    !commented out the following block of code as the subroutine does not seem to be used anywhere
   !     !    !calculate the contribution to xi_out from the matrix product of lambda and xi
   !     !    do M = 1,sys%nex
   !     !    do N = 1,sys%nex
   !     !        do I = 1,sys%nex
   !     !            do J = 1,sys%nex
   !     !                do R = 1,sys%nex
        !                     do S = 1,sys%nex
        !                         ham%xi_out(I,J,M,N) = ham%xi_out(I,J,M,N) + ham%xi(R,S,M,N)*ham%lambda(I,J,R,S)
        !                     end do
        !                 end do
        !             end do
        !    !         end do
        !    !     end do
        !    ! end do
        !!end subroutine compute_xi_out
        !subroutine compute_xi_dir_lambda_e()
        !
        !    integer :: I,J,M,N,R,S
        !    do M = 1,sys%nex
        !        do N = 1,sys%nex
        !            do I = 1,sys%nex
        !                do J = 1,sys%nex
        !                    do R = 1,sys%nex
        !                        do S = 1,sys%nex
        !                            ham%xi_dir_lambda_e(I,J,M,N) = ham%xi_dir_lambda_e(I,J,M,N) &
        !                                                          + ham%xi(R,S,M,N)*ham%lambda_e(I,J,R,S)
        !                        end do
        !                    end do
        !                end do
        !            end do
        !        end do
        !    end do
        !end subroutine compute_xi_dir_lambda_e
        !
        !subroutine compute_xi_in_lambda_e()
        !
        !    integer :: I,J,M,N,R,S
        !    do M = 1,sys%nex
        !        do N = 1,sys%nex
        !            do I = 1,sys%nex
        !                do J = 1,sys%nex
        !                    do R = 1,sys%nex
        !                        do S = 1,sys%nex
        !                           ham%xi_in_lambda_e(I,J,M,N) = ham%xi_in_lambda_e(I,J,M,N) &
        !                                                         + ham%xi_in(R,S,M,N)*ham%lambda_e(I,J,R,S)
        !                        end do
        !                    end do
        !                end do
        !            end do
        !        end do
        !    end do
        !end subroutine compute_xi_in_lambda_e
        !
        !subroutine compute_xi_out_lambda_e()
        !
        !    integer :: I,J,M,N,R,S
        !    do M = 1,sys%nex
        !        do N = 1,sys%nex
        !            do I = 1,sys%nex
        !                do J = 1,sys%nex
        !                    do R = 1,sys%nex
        !                        do S = 1,sys%nex
        !                            ham%xi_out_lambda_e(I,J,M,N) = ham%xi_out_lambda_e(I,J,M,N) &
        !                                                           + ham%xi_out(R,S,M,N)*ham%lambda_e(I,J,R,S)
        !                        end do
        !                    end do
        !                end do
        !            end do
        !        end do
        !    end do
        !
        !end subroutine compute_xi_out_lambda_e
        !
        !subroutine compute_E_MNRS_lambda_lambda_e()
        ! integer :: I,J,M,N,R,S
        !    do M = 1,sys%nex
        !        do N = 1,sys%nex
        !            do I = 1,sys%nex
        !                do J = 1,sys%nex
        !                    do R = 1,sys%nex
        !                        do S = 1,sys%nex
        !                           ham%E_MNRS_lambda_lambda_e(I,J,M,N) = ham%E_MNRS_lambda_lambda_e(I,J,M,N) &
        !                                + ((exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N) + exciton_sys%eigenvalues(R) + exciton_sys%eigenvalues(S))&
        !                                * ham%lambda(R,S,M,N)*ham%lambda_e(I,J,R,S))
        !                        end do
        !                    end do
        !                end do
        !            end do
        !        end do
        !    end do
        !
        !end subroutine compute_E_MNRS_lambda_lambda_e
        !
        !
        !subroutine compute_lambda_lambda_e()
        !
        !     integer :: I,J,M,N,R,S
        !
        !     do M = 1,sys%nex
        !        do N = 1,sys%nex
        !            do I = 1,sys%nex
        !                do J = 1,sys%nex
        !                    do R = 1,sys%nex
        !                        do S = 1,sys%nex
        !                            ham%lambda_lambda_e(I,J,M,N) = ham%lambda_lambda_e(I,J,M,N) &
        !                                                          + ham%lambda(R,S,M,N)*ham%lambda_e(I,J,R,S)
        !                        end do
        !                    end do
        !                end do
        !            end do
        !        end do
        !    end do
        !
        !
        !end subroutine compute_lambda_lambda_e
        
end module compute_lambda_and_xi
    
