	program test
	implicit none
        integer nj,ns,r,kflag
        parameter (r=12,ns=128,kflag=-1)
        parameter (nj=ns)
        integer i,iflag,xsub(nj),ier,num,mm
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(r,nj),V1(r,ns),U2(r,nj),V2(r,ns)
        real*16 x(nj),pi,time1,time2
        real*16 arr(4)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(r,nj),V(r,ns),c(ns),S(nj)
        complex*16 fk(nj)
        complex*16,allocatable :: VV(:,:)
        real*8 x1(nj),eps,error
        double complex in1, out1
        dimension in1(ns), out1(ns)
	integer*8 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001

        eps=1E-12
        num=10000
        open(unit = 10,file = 'Ur2.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr2.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui2.txt')
        read(10,*) U2
        open(unit = 10,file = 'Vi2.txt')
        read(10,*) V2
        call dfftw_plan_dft_1d(plan,nj,in1,out1,FFTW_FORWARD,0)
        print *,'start 1D type 2 testing:'
        print *,'nj  =',nj,'ns  =',ns
        print *,'rank r = ',r
        print *,'eps             =',eps
      
        U=dcmplx(U1,U2)
        V=dcmplx(V1,V2)
        V=conjg(V)
        if (kflag .lt. 0) then
           mm=floor(ns/2.0+0.6)
           allocate(VV(r,mm))
           VV=V(:,1:mm)
           V(:,1:mm)=V(:,mm+1:ns)
           V(:,mm+1:ns)=VV
        endif
           
        !print *,V(1,:)
        !print *,U(1,:)
        do i = 1,nj
           x(i) = i*pi/8
        enddo
        xsub=mod(floor(x+0.5),ns)+1
        do i = 1,ns
           c(i) = exp(-dcmplx(0,1)*i/ns)
        enddo
        !print *,c(-64:-59)
        do i = 1,nj
           x1(i) = i*pi*2*pi/(8*nj)
        enddo


        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1dIIapp(nj,plan,c,U,V,xsub,ns,kflag,r,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(1:5)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1d2f90(nj,x1,fk,-1,eps,ns,c,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
        !print *,'ier=',ier
        !print *,'fk(1:5)=',fk(1:5)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
        error=sqrt(real(sum((S-fk)*conjg(S-fk))/
     &  sum(S*conjg(S))))
        print *,' relative error= ',error
        call dfftw_destroy_plan(plan)
        

	end program

