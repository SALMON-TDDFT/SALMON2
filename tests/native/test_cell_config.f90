program probe
 use structures,only:s_cell_config,s_cell_geometry
 use initialization_sub,only:resolve_cell
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(s_cell_config)::a,b
 type(s_cell_geometry)::ga,gb,before
 integer::status
 a%lengths=[8d0,12d0,16d0];a%points=[8,12,16]
 a%vectors=0d0;a%vectors(1,1)=10;a%vectors(2,2)=20;a%vectors(3,3)=30
 call resolve_cell(a,ga,status)
 if(status/=0.or.any(ga%lengths/=a%lengths).or.any(ga%spacing/=1d0))stop 1
 b=a;b%lengths=0d0
 call resolve_cell(b,gb,status)
 if(status/=0.or.any(gb%lengths/=[10d0,20d0,30d0]))stop 2
 b%vectors=0d0;b%spacing=[2d0,3d0,4d0]
 call resolve_cell(b,gb,status)
 if(status/=0.or.any(gb%lengths/=[16d0,36d0,64d0]))stop 3
 b=a;b%vectors(:,1)=[8d0,1d0,0d0];b%vectors(:,2)=[0d0,12d0,0d0];b%vectors(:,3)=[0d0,0d0,16d0]
 call resolve_cell(b,gb,status)
 if(status/=0.or.gb%orthogonal.or.any(gb%primitive/=b%vectors))stop 4
 if(gb%lengths(1)/=sqrt(65d0).or.any(ga%lengths/=a%lengths))stop 5
 before=gb;b%spacing(1)=ieee_value(0d0,ieee_quiet_nan)
 call resolve_cell(b,gb,status)
 if(status==0.or.any(gb%primitive/=before%primitive))stop 6
 b=a;b%points=0
 call resolve_cell(b,gb,status)
 if(status==0.or.any(gb%spacing/=before%spacing))stop 7
 b=a;b%vectors(:,1)=[1d0,1d0,0d0];b%vectors(:,2)=[2d0,2d0,0d0];b%vectors(:,3)=[0d0,0d0,1d0]
 call resolve_cell(b,gb,status)
 if(status==0.or.any(gb%primitive/=before%primitive))stop 8
 print *,'Explicit cell resolution PASS'
end program
