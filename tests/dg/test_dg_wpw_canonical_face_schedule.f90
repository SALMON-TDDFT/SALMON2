program test_dg_wpw_canonical_face_schedule
  use dg_wpw_canonical_face_schedule,only:s_dg_wpw_canonical_face_record,&
    build_dg_wpw_canonical_face_schedule
  implicit none
  type(s_dg_wpw_canonical_face_record),allocatable::faces(:)
  integer::info
  call build_dg_wpw_canonical_face_schedule(1,[2,1,1],faces,info)
  if(info/=0.or.size(faces)/=2)error stop 1
  if(any([(faces(info)%neighbor_fragment/=2,info=1,size(faces))]))error stop 2
  if(.not.any([(faces(info)%side_from_local==-1,info=1,size(faces))]).or.&
     .not.any([(faces(info)%side_from_local==1,info=1,size(faces))]))error stop 3
  if(count([(faces(info)%canonical_owner,info=1,size(faces))])/=2)error stop 4
  call build_dg_wpw_canonical_face_schedule(2,[3,1,1],faces,info)
  if(info/=0.or.size(faces)/=2)error stop 5
  if(count([(faces(info)%canonical_owner,info=1,size(faces))])/=1)error stop 6
  if(.not.any([(faces(info)%neighbor_fragment==1.and..not.faces(info)%canonical_owner,&
    info=1,size(faces))]))error stop 7
  if(.not.any([(faces(info)%neighbor_fragment==3.and.faces(info)%canonical_owner,&
    info=1,size(faces))]))error stop 8
  print '(a)','PASS bounded canonical-face ownership schedule'
end program
