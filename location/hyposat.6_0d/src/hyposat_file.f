cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function file_check(name)
c
c     function checks if file exists and if eventually
c     the file must be searched in the directory
c     defined by environment variable HYPOSAT_DATA
c
c     26 October 2000
c
c     last changes: 10 September 2013
c

      character*(*) name
      character*120 path, envvar, file, file_check

      logical exst

      exst = .false.

      file_check = ' '

      if(name.eq.'_') return

      ifil = len_trim(name)
      inquire(file=name(1:ifil),exist=exst)

      if(exst) then

         file_check = name(1:ifil)

      else

c
c        Start of a compiler-specific block
c        Get name of directory path from UNIX environment variable:
c        HYPOSAT_DATA. This part has eventually be changed for
c        other machines!!!
c
         envvar = 'HYPOSAT_DATA'
         CALL get_environment_variable(ENVVAR, path)

         file = path(1:len_trim(path)) // '/' // name(1:ifil)

         ifil = len_trim(file)

         inquire(file=file(1:ifil),exist=exst)

         if(exst) then

            file_check = file(1:ifil)

         else

            write (*,'('' Major ERROR: Cannot find to open: '',a)') 
     *            name
            write (*,'('' Major ERROR: Cannot find to open: '',a)') 
     *            file

            write (*,'(''Set correct path name with environment '',
     *            ''variable HYPOSAT_DATA !'')')
            stop

         endif

      ENDIF

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function file_checkpara()
c
c     function checks if a file called 'hyposat-parameter' either 
c     exists in the current directory, or if this file exists in the 
c     directory defined by the environment variable HYPOSAT_PARA
c     (i.e., HYPOSAT_PARA has only the "PATH"), or if the file
c     exists defined by HYPOSAT_PARA (i.e., HYPOSAT_PARA contains
c     the whole file name).
c
c     The search order is:
c
c        1) hyposat-file in the current directory
c        2) hyposat-file in the directory defined by HYPOSAT_PARA
c        3) HYPOSAT_PARA contains the whole file name
c
c     5 November 2003
c
c     last changes: 10 September 2013
c

      character*120 path, envvar, file_checkpara, file1, file2

      logical exst

      exst = .false.

      file_checkpara = ' '

      inquire(file='hyposat-parameter',exist=exst)

      if(exst) then

         file_checkpara = 'hyposat-parameter'

      else

c
c        Start of a compiler-specific block
c        Get name of directory path from UNIX environment variable:
c        HYPOSAT_PARA. This part has eventually be changed for
c        other machines!!!
c
         envvar = 'HYPOSAT_PARA'
         CALL get_environment_variable(ENVVAR, path)

         file1 = trim(path) // '/hyposat-parameter'

         inquire(file=trim(file1),exist=exst)

         if(exst) then

            file_checkpara = trim(file1)

         else

            file2 = trim(path)

            inquire(file=trim(file2),exist=exst)

            if(exst) then

               file_checkpara = trim(file2)

            else

               write (*,'('' Major ERROR: Cannot find to open: '',
     *            ''./hyposat-parameter'')')
               write (*,'('' Major ERROR: Cannot find to open: '',a)') 
     *            file1
               write (*,'('' Major ERROR: Cannot find to open: '',a)') 
     *            file2
               write (*,'(''Set correct path/file-name with '',
     *            ''environment variable HYPOSAT_PARA !'')')

            endif
         endif

      ENDIF

      return
      end
