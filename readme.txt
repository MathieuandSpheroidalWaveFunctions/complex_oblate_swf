                               complex_oblate_swf

  Coblfcn is available as both a subroutine version provided as the
  module complex_oblate_swf and a stand alone version coblfcn. It was
  originally developed by arnie lee van buren about 6 years years ago
  and has been updated and improved several times since then.

  Table of Contents
  1. Purpose
  2. Introduction
  3. Input and Output
  4. Accuracy of results
  5. Obtaining the expansion d coefficients
  6. Obtaining the eigenvalues


 1. Purpose

  To calculate the first and second kind oblate radial functions r1
  and r2 and their first derivatives r1d and r2d for a given order m,
  a range of degrees l beginning at l = m and for a specific complex size
  parameter c and shape parameter x. To calculate the first kind oblate
  angular functions and their first derivatives with respect to eta for
  the same values of m, l and c and for a set of values of the angular
  coordinate eta. The subroutine version complex_oblate_swf calculates
  values for a single input value of m. The stand alone version calculates
  values for a range of values of m.

 2. Introduction

  Coblfcn is written in free format fortran. It is designed around the
  maximum number of decimal digits ndec and the maximum exponent nex
  available in real arithmetic. Procedures used in oblfcn allow for
  exponents much larger than nex since the resulting floating point
  radial function values are given as a characteristic and an integer
  exponent.

  Coblfcn can be run in double precision, quadruple precision or a hybrid
  where the Bouwkamp procedure to refine the eigenvalues is run in quadruple
  precision while the remainder of the calculations are performed in double
  precision. In the latter case, coblfcn switches to quadruple precision for
  the Bouwkamp procedure whenever double precision fails to provide highly
  accurate eigenvalues. See the discussion below about when to choose which
  arithmetic. The choice is set in the module param provided in the github
  repository. If this is not available, then create param as follows:

    module param
    integer, parameter :: knd = selected_real_kind(8)
    integer, parameter :: knd1 = selected_real_kind(8) 
    logical, parameter :: debug = .true.
    logical, parameter :: warn = .true.
    logical, parameter :: output = .true.
    logical, parameter :: suffix = .true.
    end module param

  Set the value of knd in the parenthesis to either 8 for double precision
  or 16 for quadruple precision arithmetic. Set the value of knd1 to be used
  for the Bouwkamp procedure. Note that if knd = 16, knd1 should also be 16.
  The logicals in param are described below in the discussion of the output
  files.

  Some computers may have more than 8 bytes for double precision
  data and more than 16 bytes for quadruple precision data or may use
  kind values that do not correspond to the number of bytes. In this
  case just use the appropriate integers for the kind parameters in
  module param. Also change the values of kindd and kindq set in
  statement 5 below below the comments section to the kind values for
  double precision data and quadruple precision data, respectively.

  A description of the methods used in coblfcn is provided in the manuscript
  'Calculation of oblate spheroidal wave functions with complex argument,'
  available at arXiv.org, identifier 2009.01618, August 2020. Coblfcn was
  developed around oblfcn which obtains spheroidal function values for real
  values of the size parameter. A description of the methods used in
  oblfcn is provided in the article 'Accurate calculation of oblate
  spheroidal wave functions,' available at arXiv.org, identifier
  1708.07929, August 2017, revised September 2019.

  Coblfcn provides function values for c complex = real(c) + i aimag(c)
  = cr + i ci, where the imaginary part ci often accounts for losses in
  wave propagation. Ci is assumed positive in coblfcn. If the user has
  a negative value for ci, just run coblfcn with ci positive instead
  and take the complex conjugate of the results, including the function
  values, eigenvalues, expansion coefficients, and normalization
  factors.

  3. Input and Output

  Following is a description of the input and output parameters in the
  call statement for the subroutine version. After that will be a
  description of the the input and output files associated with the
  stand alone version. Note that these output files, if desired, can
  also be obtained when running the subroutine version. See comments about
  this below.

  A sample input and resulting output from oblfcn is provided by the
  files coblfcndat (text version of the input file coblfcn.dat for the
  stand alone version), coblfort20 (text version of the output file
  fort.20 giving the resulting radial functions) and coblfort30 (text
  version of the output file fort.30 giving the resulting angular
  functions).

  Subroutine Version of coblfcn

    subroutine oblfcn(c,m,lnum,ioprad,x,iopang,iopnorm,narg,arg, &
                      r1c,ir1e,r1dc,ir1de,r2c,ir2e,r2dc,ir2de,naccr, &
                      s1c,is1e,s1dc,is1de,naccs)

        complex(knd), intent(in)   ::  c
        real(knd), intent (in)     ::  x, arg(narg)
        integer, intent (in)       ::  m, lnum, ioprad, iopang, iopnorm, narg
        complex(knd), intent (out) ::  r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                                       s1c(lnum, narg), s1dc(lnum, narg)
        integer, intent (out)      ::  ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                                       is1e(lnum, narg), is1de(lnum, narg), & 
                                       naccr(lnum), naccs(lnum, narg)

      Input and output parameters appearing in the subroutine call
      statement are defined below:

          c      : desired complex value of the size parameter (= kd/2,
                   where k =  the complex wavenumber and d = interfocal
                   length) [real(knd)]
          m      : desired value for the order m (integer)
          lnum   : number of values desired for the degree l equal
                   to m, m + 1, m + 2, ..., m + lnum - 1 (integer)
                   if lnum is less than 2*[real(c)+aimag(c)]/pi it
                   should be an even integer.
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed
          x      : value of the radial coordinate x [real(knd)]
          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed
          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner and
                      Schafke normalization scheme.] This norm
                      becomes very large as m becomes large. The
                      angular functions are computed below as
                      a characteristic and an exponent to avoid
                      overflow.
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.
                      This is very useful since it removes the
                      need to calculate a normalization factor
                      when using the angular function values given
                      here. It also eliminates any chance for
                      overflow when the characteristics and exponents
                      are combined to form the angular functions.
          narg   : number of values of the angular coordinate eta for
                   which angular functions are calculated (integer)
          arg    : vector containing the narg values of eta for which
                   angular functions are desired [real(knd)]
          r1c/   : vectors of length lnum containing the complex characteristics
          r1dc     for the radial functions of the first kind r1 and their
                   first derivatives [complex(knd)]
          ir1e/  : integer vectors of length lnum containing the
          ir1de    exponents corresponding to r1c and r1dc
          r2c/   : vectors of length lnum containing the complex characteristics
          r2dc     for the radial functions of the second kind r2 and their
                   first derivatives [complex(knd)]
          ir2e/  : integer vectors of length lnum containing the
          ir2de    exponents corresponding to r2c and r2dc
          naccr  : integer vector of length lnum containing the estimated
                   accuracy of the radial functions
          s1c,   : two-dimensional arrays s1c(lnum,narg) and
          s1dc     s1dc(lnum,narg) that contain narg calculated
                   compelx characteristics for the angular functions and
                   their first derivatives for each of the lnum
                   values of l [complex(knd)]
                   For example, s1(10,1) is the characteristic
                   of the angular function for l = m +10 -1 and
                   for the first value of eta given by arg(1)
          is1e,    integer arrays is1e(lnum,narg) and is1de(lnum,narg)
          is1de    containing the exponents corresponding to s1c and
                   s1dc
          naccs  : two-dimensional array naccs(lnum,narg) that contains
                   narg estimated accuracy values for the angular functions
                   for each of the lnum values of l

  Stand alone Version of oblfcn

     Input Data

     Input parameters are read from unit 1 in the file coblfcn.dat
     assumed to be in the directory of coblfcn.f90. coblfcn.dat
     contains the following lines of data:

       line 1:
          mmin   : minimum value for m. (integer)
          minc   : increment for m. (integer)
          mnum   : number of values of m. (integer)
          lnum   : number of values of l [l=m, l=m+1,
                   ..., l=m+lnum-1]. (integer)
                   if lnum is less than 2*[real(c)+aimag(c)]/pi
                   it should be an even integer. If lnum is chosen
                   to be odd, oblfcn will increase lnum by one.

       line 2:
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed

          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed

          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner-
                      Schafke normalization scheme.]
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.

       line 3:
          c      : value of the complex size parameter (= kd/2, where k =
                   the complex wavenumber and d = interfocal length) [complex(knd)]
          x      : value of the radial coordinate x [real(knd)]
                   (any value can be entered if ioprad = 0)

       line 4:
          ioparg : (integer)
                 : =0 if both arg1 and darg are angles in degrees
                 : =1 if arg1 and darg are dimensionless values of eta

          arg1   : first value for the angle coordinate (in degrees
                   or dimensionless if eta) for which angular
                   functions are to be computed. [real(knd)]

          darg   : increment used to calculate additional desired
                   arguments for angular functions. [real(knd)]

          narg   : number of desired angle arguments. (integer)
                   (line 4 is not read when iopang = 0)

     Output files

     These output files are also available using the subroutine version
     complex_oblate_swf. Generation of each of the files is controlled by a
     logical specified in the module param. False suppresses the output file
     and true enables it. The logical debug controls fort.30 and fort.40,
     the logical output controls fort.20 and fort.30 and warn controls fort.60.
     The logical suffix controls whether the accuracy estimates given in fort.20
     are followed by a letter designating how the accuracy was determined. 'w'
     indicates it is based on the Wronskian and 'e' indicates it is based on
     subtraction errors involved in the calculations. Setting suffix = false
     suppresses the letter. 

   fort.20

     This file contains values for all radial functions that have
     been calculated.
     The first line in the file contains the values for x, c, and
     m, formatted as follows (see statements 260 and 265 in
     subroutine main):

                x      : e23.14 in 64 bit arithmetic; e38.30
                       : in 128 bit arithmetic
                c      : (e23.14,e23.14) in 64 bit arithmetic;
                       : (e38.30,e38,30) in 128 bit arithmetic
                m      : i5

     Each subsequent line in fort.20 contains radial functions
     for given values of l. The first line contains values for l = m,
     the next for l=m+1 and continuing to l=m+lnum-1. The radial
     functions are preceeded by the value of l and followed by the
     accuracy, equal to the estimated number of accurate decimal digits
     in the radial functions as measured using the Wronskian or
     estimated using subtraction errors in the calculations. [see
     comments below regarding naccr]

       The output and corresponding format for each line is as follows
       (see statements 1340, 1350, and 1380 in main).

         l            : value for l (i5)
         r1c(l-m+1)   : complex characteristic of the oblate radial
                        function of the first kind r1 for
                        the given value of l (f17.14,f17.14)
         ir1e(l-m+1)  : exponent of r1 (i6)
         r1dc(l-m+1)  : complex characteristic of the first derivative
                        of the oblate radial function of the
                        first kind r1d (f17.14,f17.14)
         ir1de(l-m+1) : exponent of r1d (i6)
         r2c(l-m+1)   : complex characteristic of the oblate radial
                        function of the second kind r2 for
                        the given value of l (f17.14,f14.7)
         ir2e(l-m+1)  : exponent of the oblate radial function of
                        second kind (i6). If the exponent for any
                        of the functions is greater than 9999,
                        the format can be increased to i7 or
                        higher. Note that the procedures used in
                        this program allow for exponents much
                        larger than those allowed in the arithmetic
                        used for the calculations since the
                        floating point function values provided
                        are given as a characteristic and an
                        integer exponent. Use of ratios in
                        calculating the functions eliminates
                        overflow and underflow during the
                        calculations.
         r2dc(l-m+1)  : complex characteristic of the first derivative
                        of the oblate radial function of second kind
                        r2d (f17.14,f17.4)
         ir2de(l-m+1) : exponent of r2d (i6). [See comment above
                        for ir2e.]
         naccr(l-m+1) : estimated accuracy: usually equal to the
                        number of decimal digits of agreement
                        between the theoretical Wronskian and the
                        calculated Wronskian (i2). This is
                        indicated by the letter w following the
                        integer naccr. r1 and r1d are usually very
                        accurate so that naccr relects the accuracy
                        of r2 and r2d.

                        Several situations are described below where
                        the Wronskian can not be used to estimate the
                        accuracy. Other factors are used here to
                        estimate naccr. This is indicated by using
                        the letter e instead of w following the value
                        for naccr in fort.20.

                        Sometimes near the breakpoint value for l,
                        the Legendre function expansions for r2 and r2d
                        converge with acceptable subtraction error but
                        the leading coefficient for the series with the
                        P Legendre function is highly inaccurate. Here
                        the Wronskian can often be used for small ci to
                        obtain improved accuracy for the coefficient.
                        The accuracy is estimated using subtraction
                        errors involved in the calculations and the
                        estimated accuracy of the eigenvalue in those
                        cases where the Bouwkamp eigenvalue routine
                        does not converge fully.

                        When ci becomes large, there are often values
                        of l where both r1 and r2 and their first
                        derivatives are large in magnitude. Here values
                        for r2 and r2d are given by -i times values of
                        r1 and r1d. The accuracy is estimated using
                        the magnitudes of r1 and r1d, the magnitude
                        of the theoretical value of the wronskian, and
                        the estimated accuracy of r1 and r1d.

                        When x is equal to zero, the Wronskian is used
                        to obtain r2d when l-m is even and r2 when l-m
                        is odd. Here the accuracy naccr is approximated
                        using the estimated accuracy of the eigenvalue
                        and the estimated accuracy of the joining
                        factor used in using the Legendre function
                        expansion to compute r2 and r2d. The accuracy
                        estimate for the other functions is given by
                        the estimated accuracy of the eigenvalue minus
                        one digit. When r2 or r2d are less accurate due
                        to large subtraction errors in the calculation
                        of the joining factor, they are often somewhat
                        small in magnitude, especially for small values
                        of ci. Here the reduced accuracy is not likely
                        to result in lower accuracy in solutions to
                        problems utilizing these functions. When the
                        estimated accuracy is zero, the functions are
                        set equal to zero.

                        There are two other situations where the value
                        for naccr can be an estimated value: (1) When
                        the function values for r2 and r2d are
                        obtained using the integral expressions for
                        small values of x, the accuracy given by the
                        Wronskian value is adjusted downward to account
                        for the likelihood that the accuracy for r2 for
                        l - m even and r2d for l - m odd is lower than
                        the Wronskian estimate. (2) The accuracy using
                        the traditional Legendre expansion is also
                        estimated and the accuracy estimate provided
                        by coblfcn is taken to be the smaller of this
                        and the Wronskian estimate. In both cases, the
                        accuracy is designated by w rather than e since
                        the Wronskian accuracy sets an upper bound.

   fort.30

     This file fort.30 contains values for all angular functions
     that have been calculated. Its first line contains the values
     for c and m, formatted as follows (see statements 50 and 55 in
     subroutine main).

                c      : (e23.14,e23.14) in 64 bit arithmetic;
                       : (e38.30,e38,30) in 128 bit arithmetic
                m      : i5

     The second line in fort.30 contains the value for the first l
     (=m), formatted as follows (see statement 270 in subroutine
     main):

                l      : i6

     This is followed by a series of narg lines. Each line contains
     a desired value of angle (ioparg = 0) or angular coordinate eta
     (ioparg =1) followed by the corresponding angular functions and
     accuracy. Specific output and format for each line is as follows.

        for iopang = 1:

               arg    : for ioparg = 0, angle in degrees (f17.14; see
                        statement 1440 in subroutine main)
            or barg   ; for ioparg = 1, angular coordinate eta
                        (f17.14; see statement 1460 in subroutine main)
               s1c    : complex characteristic of the oblate angular
                        function of first kind (f17.14,f17.14); see
                        statement 1460 in subroutine main)
               is1e   : exponent of the oblate angular function of
                        first kind (i5; see statement 1460 in
                        subroutine main)

        for iopang = 2, each line also includes:
               s1dc   : complex characteristic of the first derivative
                        of the oblate angular function of first kind
                        (f17.14,f17.14); see statement 1470 in
                        subroutine main)
               is1de  : exponent of the first derivative of the
                        oblate angular function of first kind (i5;
                        see statement 1470 in subroutine main)

        for iopang = 1:
               naccs  : accuracy: estimate of the number of decimal
                        digits of accuracy in the angular function
                        (i2). It is a conservative estimate based on
                        the calculated subtraction error in the
                        Legendre function series for the angular
                        functions and the estimated accuracy of the
                        normaliation factor. When the accuracy estimate
                        is equal to 0, the corresponding angular
                        functions are set equal to zero. (i2; see
                        statements 1460 and 1470 in subroutine main).
                        The calculated angular functions tend to be
                        less accurate the larger the value of cr, the
                        smaller the value of l - m (for values less
                        than approximately 2c divided by pi), and the
                        closer eta is to zero (i.e., the closer theta
                        is to 90 degrees).

        for iopang = 2:
              naccds  : accuracy: includes an estimate of the number
                        of decimal digits of accuracy in the first
                        derivative of the angular function (i2)   

   fort.40 and fort.50

     These files are diagnostic files that contain information
     about specific techniques used and numbers of terms required
     for the radial function and angular function calculations,
     respectively. They are annotated and should be self
     explanatory.

   fort.60

     This file may be of interest to the user, especially when using
     this program outside the range for which coblfcn is expected to
     provide useful results. Whenever the estimated accuracy falls below
     a designated integer value during the running of this program,
     the associated values of x, c, m, and l are written to fort.60.
     The integer is currently set equal to 6 in the write statements
     for this file found after the line numbered 1405 in subroutine main.
     Whenever coblfcn determines that two or more of the eigenvalues of
     the same parity are duplicates, indicating a failure of the Bouwkamp
     routine to converge to the correct eigenvalue, the values of m,l and
     c where this occurs are written to fort.60.

  4. Accuracy of Results

  The following discussion is provided to help the user choose which
  arithmetic option to use. If the compiler does not support quadruple
  precision arithmetic, then the only option is to use double precision.
  If the compiler does support quadruple precision arithmetic, then it
  is recommended that the hybrid version be used instead of the entirely
  double precision version whenever double precision is expected to
  provide adequate accuracy. Use of quadruple precision arithmetic
  for the Bouwkamp procedure increases the run time somewhat but the
  program is still reasonably fast and the results are more accurate.
  The improvement in accuracy increases as ci increases. If the input
  parameters are outside those for which the hybrid version is expected
  to provide useful results, then the only choice is to use quadruple
  precision if it is available. Note that this will increase the run
  time, although coblfcn runs reasonably fast in quadruple precision
  when cr is not very large.

  Coblfcn was tested extensively using a laptop pc and a Fortran
  compiler that provides approximately 15 decimal digits in double
  precision arithmetic and approximately 31 digits in quadruple
  precision arithmetic. If the user's computer provides a different
  number of digits for either double or quadruple precision, the
  following estimates might need to be adjusted up or down depending
  on whether more or fewer digits are provided.

  The estimated accuracy of the resulting function values is given
  below in terms of decimal digits. It is often obtained using the
  comparison of the theoretical value for the Wronskian and its value
  obtained using the computed radial function values. But it can be an
  estimate based on subtraction errors in the calculation. See the
  discussion for the accuracy integer naccr in the section for fort.20.

  Testing included values of x ranging from 0.000001 to 10 as well as
  the special case of x = 0, values for cr up to 5000, and values of ci
  up to 200. Testing for both the double precision and hybrid versions
  included all values of the order m from 0 to 200 and from 210 to 1000
  in steps of 10. Testing for the quadruple precision version was the
  same for cr less than about 1000. Testing for larger cr included
  values of m from 0 to 200 in steps of 10 and from 250 to 1000 in
  steps of 50. For all versions the values of the degree l ranged from
  m to m + lnum -1, where lnum was chosen sufficiently large that the
  magnitudes of r1 and r1d were less than 10**(-300).

  An integer called minacc is used in coblfcn. This designates the
  minimum number of accurate digits desired for the radial spheroidal
  functions of the second kind. The value of minacc controls which
  methods are used to calculate the radial functions. Minacc is set
  equal to 8 for double precision. It is recommended that this not
  be changed. For quadruple precision minacc is set equal to 15 digits
  for values of ci up to 20. This should provide 15 or more digits of
  accuracy here. For ci > 20, minacc is set equal to 8 digits. Minacc
  can be increased in this case but higher accuracy might not always
  be achieved. Also the computation time will likely go up with larger
  values of minacc. Minacc is set below following these introductory
  comment statements.

  In the following discussion the term useful results means that the
  estimated accuracy for the radial functions observed during testing
  never fell below 5 decimal digits unless otherwise stated. I expect
  there are many applications where occasional 5 digit results are
  acceptable. Possibly even an isolated 4 digit result is acceptable.
  Note that there is no guarantee that the estimated accuracy for the
  values for parameter values other than those that I tested will be
  as high as I report below. The discussion below will focus on x
  unequal to zero. It is expected that function values for x = 0 will
  be useful for the same values of c and m that useful results are
  obtained for at least one value of x. Here both the radial function
  of the first kind r1 and its first derivative will be very accurate.
  Whenever values of the radial function of the second kind r2 for
  l - m even are less accurate than 5 digits, they are expected to be
  proportionally smaller in magnitude than r1. Similarly for the first
  derivatives of r2 and r1 when l - m is odd.

  Double precision arithmetic

  Using double precision arithmetic, coblfcn provides useful results
  for ci (the imaginary part of c) as large as 10, for cr (the real
  real part of c) as large as 5000, for m up to at least 1000 and for
  all tested values of x down to 0.000001. When ci is less than about
  5, coblfcn provides results with an accuracy similar to that
  provided by the program oblfcn for real c. Extensive testing for
  ci = 10 showed that the radial functions of the first kind r1 and
  its first derivative r1d are almost always accurate to 10 or more
  decimal digits. For all values of x except zero the radial functions
  r2 and r2d are usually accurate to 8 or more decimal digits but
  accuracies lower than this were seen, especially for larger cr and
  smaller x. Nearly all accuracies less than 8 digits occurred near but
  somewhat below the so-called breakpoint. [The breakpoint is defined
  as the minimum value of l such that the magnitudes of r2 and r2d
  begin to increase with increasing l while the magnitudes of r1 and
  r1d begin to decrease. An approximate value for the breakpoint
  for small values of m is given by l = 2*(cr+ci)/pi]. No 5 digit
  results were seen for cr up to 200 or so. Only a few 5 digit results
  were seen, even for cr = 5000. See the discussions below about the
  estimated accuracy of the angular functions.

  Similar testing when ci was increased to 12 showed a few more 5 digit
  results. At least 5 digits were obtained for cr up to 5000 for all of
  the values of x down to 0.000001.

  Testing with ci = 15 showed at least 5 digits of accuracy for cr up to
  2000 for all values of x down to 0.000001. Testing at cr = 5000 showed
  poor results with accuracies as low as 0 digits for m = 0. This was due
  to a poor matrix starting value for the eigenvalue together with a lack
  of convergence of the Bouwkamp procedure, even when using quadruple
  precision for this procedure. 

  When ci = 20, there were 5 or more digits of accuracy for all tested
  values of x when cr <= 70. There was one estimated accuracy of 4 digits
  for x = .01, but it was found to actually be accurate to 5 digits when
  compared to quadruple precision results. There were 5 or more digits
  of accuracy at x >= 0.1 for cr = 150 and at x >= 0.2 for cr up
  to 2000. Testing showed duplicated eigenvalues of the same parity for
  some values of m when cr = 5000. See the discussion of fort.60 above.

  When ci = 25, there were at least 5 digits of accuracy for x >= 0.2
  when cr <= 100, for x >= 0.3 when cr = 150 and for x >= 0.4 when cr
  = 200.

  Testing for yet higher values of ci showed a continued increase in
  the minimum value of x and a decrease in the maximum value of cr for
  which useful function values were obtained. If the user is interested
  in values of x somewhat larger than 0.3 together with moderate to
  small values of cr, then the double precision, or more likely the
  hybrid version, of coblfcn may be useful when ci is greater than 25.
  Otherwise, it will be necessary to use the quadruple precision
  version. Note that coblfcn runs reasonably fast in quadruple precision
  for small to moderate values of cr. If you use double precision here,
  it is recommended that the file fort.60 described below be used to
  assure that you are obtaining the accuracy you need and that there are
  no repeated eigenvalues of the same parity. Testing showed the appearance
  of duplicated eigenvalues for ci = 45 for some values of m when cr was
  only 200.

  Quadruple precision arithmetic

  When lower accuracy occurs using double precision arithmetic, much
  higher accuracy can be obtained using quadruple precision, but with
  a considerable increase in computation time, especially for very large
  cr.

  Testing for ci = 30 showed that useful results were obtained for all
  values of x down to 0.000001 for cr up to at least 5000. There were
  a few isolated 5 digit results. See the discussions of estimated
  accuracies given below for the special case of x = 0 and for the angular
  functions.

  Testing for ci = 40 showed that useful results were obtained for all
  values of x down to 0.000001 for cr up to about 4000. There were some
  scattered 5 digit results. Results for cr = 5000 showed a few isolated
  4 digit results occurring primarily for x from about 0.009 to 0.025.
  
  Testing for ci = 50 showed useful results for x down to 0.000001 and cr
  up to at least 600. Accuracies were usually 6 or more digits with a
  rare 5 digit result for m less than 200. Above m = 200, there were a
  number of 5 digit results and even an occasional 4 digit result for m
  between 200 and 350. There were also useful results for cr up to 2000
  when x <= 0.01 and m was no larger than about 600.

  Testing for ci = 60 showed useful results for all m for cr up to 150
  when x >= 0.01, for cr up to 300 when x >= 0.1 and for cr up to 1500
  when x >= 0.5.

  Testing for ci = 70 showed useful results for all m for cr up to 75
  when x >= 0.1, for cr up to 200 when x >= 0.2 and for cr up to 500
  when x >= 0.5, although a rare 4 digit result occurred sometimes.

  Testing for ci = 80 showed useful results for all m for cr up to 30
  when x >= 0.2, for cr up to 180 when x >= 0.3 and for cr up to 250
  when x >= 0.5, although a rare 4 digit result occurred sometimes.

  Testing for yet higher values of ci showed a continued increase in
  the minimum value of x and a decrease in the maximum value of cr
  for which useful function values were obtained. It is recommended
  that the user, especially for ci greater than 80, use the file
  fort.60 described below to assure that the desired accuracy is being
  obtained and that there are no repeated eigenvalues of the same
  parity.

  Sometimes either r2 or r2d will have lower accuracy than the
  estimate provided by coblfcn. This can happen for large values of
  cr, values of x less than about 0.1 and a few values of l - m
  somewhat below or near the so-called breakpoint. Here r2 for l - m
  even or r2d for l - m odd may be less accurate. However, these
  function values are usually somewhat smaller in magnitude than the
  corresponding values for r1 or r1d and make a reduced contribution
  to the numereical solution of physical problems utilizing oblate
  spheroidal functions.

     Estimated angular function accuracy

  For both choices of arithmetic, the angular functions and their first
  derivatives can lose accuracy for low l, high c, and eta near zero
  due to subtraction errors in their series calculation. However,
  their magnitude in this case is corresponding smaller than angular
  functions for higher values of l and/or eta not near zero. The loss
  in accuracy due to these subtraction errors should not adversely
  effect numerical results for physical problems using these functions.

  A second source of inaccuracy in the angular functions arises from
  subtraction errors that can occur in the calculation of their
  Meixner-Schafke normalization factor dmsnorm. Note that the loss in
  accuracy here is not in addition to other losses in accuracy for the
  angular functions but rather sets an upper limit to their accuracy.
  These errors are largest for m = 0. They occur for values of l - m
  somewhat less than the breakpoint and grow with increasing ci and to
  some extent with increasing cr. They are zero for small ci and can
  become as large as 6 digits for ci = 20 and 12 digits for ci = 50 as
  cr increases to 5000. This loss in accuracy is not likely a problem
  when using double precision arithmetic with 15 decimal digits since
  as ci becomes larger than 20, the values of cr for which the radial
  functions are accurate to 5 or more digits are progressively smaller,
  being 100 for ci = 25 and less than this for higher ci. Using double
  precision, the Meixner and Schafke normalization should be accurate
  to at least 5 digits wherever all of the radial functions for a
  given value of m are also accurate to at least 5 digits. The file
  fort.60 mentioned above can also alert the user whenever the
  estimated accuracy for the Meixner and Schafke normalization is less
  than the same integer number of decimal digits selected for alerts
  about the accuracy of the radial functions.

  For higher values of ci when using quadruple precision, the loss of
  accuracy in the normalization factor is even greater. For ci = 60,
  the loss of accuracy can be as large as 24 digits for cr = 2000 and
  25 digits for cr = 5000. For ci = 70 it is 26 digits for cr = 1000.
  And for ci = 80 it is 25 digits for cr = 400. This should not be a
  problem using 33 decimal digits since it still allows for accuracies
  of at least 5 digits for the angular functions everywhere the radial
  functions also have an accuracy of 5 or more digits.

  A third source of inaccuracy in the angular functions arises from the
  potential loss of accuracy in the eigenvalues at values of l near and
  somewhat below the breakpoint when ci is not very small and ci is
  moderate to large. This ia most likely to occur when using double
  precision for the calculations including for the Bouwkamp procedure.
  Use of quadruple precision, if available, for the Bouwkamp procedure
  will help considerably here.  

  5. Obtaining the d expansion coefficients

  The user may desire values for the d coefficients that appear in
  the expression for the angular functions as well as in many of the
  expressions used to calculate the radial functions. Ratios of
  successive d coefficients are stored in the vector enr where enr(k)
  = d(subscript 2k+ix) divided by d(subscript 2k-2+ix). The vector enr
  is calculated in the subroutine dnorm in statement 20 and passed to
  subroutine main. The number lim2 of d coefficients calculated for a
  given l is chosen to be sufficient to compute radial and angulalar
  functions for that l. The number of d coefficients needed to compute
  r1 and r1d and s1 and s1d can range from less than 100 to somewhat
  more than l/2 for large values of l. The number of d coefficients
  required to compute r2 and r2d are comparable to this unless they are
  computed using one of the Neumann function expansions. Here one can
  require up to 200000 or so coefficients when x is near 0.01. Note
  that the vector enr returned by the subroutine conver contains scaled
  ratios where the scaling has been chosen to produce a symmetric
  matrix for computing eigenvalues. The scaling factors are removed in
  subroutine dnorm to obtain the desired d coefficient ratios.

  The d coefficients themselves can be obtained starting with the value
  for d with the subscript l - m. If iopnorm is set = 0, coblfcn uses
  the Meixner_Schafke normalization scheme for the angular functions.
  Here the angular functions have the same norm as the corresponding
  associated Legendre functions. When c is real the calculation of the
  normalizing factor is accurate with no associated subraction errors.
  As ci increases and to a some extent as cr increases, the subtraction
  errors increase for some values of l. They are near zero for ci up to
  10, then increase to 6 digits at ci = 20 and to digits at ci = 50.
  See the discussion of this in the section 'Estimated angular function
  accuracy' given above. The subroutine dnorm computes d(subscript l-m)
  for this normalization and returns it to subroutine main as a
  characteristic dmlms and an exponent idmlmse. Use of an exponent
  avoids possible overflow of d(subscript l-m) for extremely large c
  and m. When the user sets iopnorm = 1 so that the angular functions
  have unit norm, the corresponding characteristic and exponent for
  d(subscript l-m) are calculated in subroutine s1 and returned to
  subroutine main as dmlms1 and idmlms1e. Corresponding values for the
  characteristic and exponent of d(subscript l-m) for the  Morse-
  Feshbach and Flammer normalizations are computed in dnorm and
  returned to main as dmlmf, idmlmfe and dmlf, idmlfe. Note that for
  c complex, all three of these normalization calculations suffer
  subtraction errors for lower values of l-m and non-small c. The
  values for d(subscript l-m) will have reduced accuracy in this
  case. When c is real only the Flammer normalization suffers
  subtraction errors in its calculation.
 
  6. Obtaining the eigenvalues

  The eigenvalues for the oblate functions are computed in subroutine
  conver and returned to main where they are stored in the vector eig(l+1).
  There is such a vector created for each value of m.
