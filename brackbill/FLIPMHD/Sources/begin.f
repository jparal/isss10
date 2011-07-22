      subroutine begin(name)
!
      use vast_kind_param, ONLY:  double
      use corgan_com_M
      implicit none
!

      character*80 :: name
!l    krd    ...device number for read input cards
!l    kpt    ...device number for print
!l    kpr    ...device number for print (switched)
!l    kpl    ...device number for storing data on disk
!l    kfm    ...device number for film/microfiche
!l    ktpin  ...device number for binary tape read
!l    ktpout ...device number for binary tape write
         krd   = 5
         kpt   = 49 
         kpl   = 15
         kfm   = 59
         ktpin =  7
         ktpout=  7
!
        open (krd, file='input')
        open (kfm, file='Initial.dat')
        open(unit=41,file='Timing.dat',status='replace')
!       open (kpt, file='Cel.out')
!
!      get the title
       read(krd,1111) name
 1111  format(a80)
!
!     get date and time:
       return
       end
