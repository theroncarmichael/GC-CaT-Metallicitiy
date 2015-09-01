;------------------------------------------------------------------------------

;Rewritten version of gcppxf.pro originally by Chris Usher, 2011
;Edited to fit DEIMOS spectra (observed by Forbes, Foster, Pastorello, 2013-09-29) by Theron Carmichael, Jun 2015

;To save velocity measurements to a file, add the necessary file writing lines to this code or somewhere in the main ppxf.pro code

pro deimos_ppxf

  ;wavelengthMask = [[8428, 8438], [8450, 8455], [8463, 8468], [8492, 8496], [8502, 8507], [8536, 8541], [8546, 8551], [8586, 8591], [8595,8599], [8610, 8687], [8757, 8781], [8823, 8839], [8865,8870]]
  wavelengthMask = [[8413, 8418], [8428, 8438], [8450, 8455], [8463, 8468], [8492, 8496], [8502, 8507], [8536, 8541], [8546, 8551], [8586, 8591], [8595,8599], [8610, 8687], [8757, 8781], [8823, 8839], [8848, 8852], [8865,8870], [8882, 8888], [8900, 8906], [8917, 8923], [8941, 8946], [8955, 8961]]

  ;The number of times each spectrum is resampled. Doing this adds artifical noise and attempts to re-fit a spectrum, resulting in
  ;better reducde chi-squared values: see line 111
  resamplings = 3

  ;Template spectra library; file_search() deals with a list of files, hence the '*.fits'
  templatefiles = file_search('SAGES/NGC1023/template/l_t*.fits', COUNT = ntemplates)

  ;Object files containing spectra
  spectrafiles = file_search('SAGES/NGC1023/spectraR/CaTrip/n1023*.fits', COUNT = sfiles)

  ;Noise files corresponding to spectra files
  noisefiles = file_search('SAGES/NGC1023/varianceR/n1023*.fits', COUNT = nofiles)

  ;Velocity list that have been fitted by pPXF once, ready to resample
  ;For first time use of pPXF, run initial guess for velocity through a range of velocity values (for loop) to find the best fit
  ;readcol, 'SAGES/NGC1023/spectraR/CaTrip/gcppxf_v_simp_1023.dat', name, veloc , disp, chi2, format = 'a,f,f,f'



  ;Run this loop for each of the spectra files in the spectrafiles library (55 in this example)
  ;Can remove for loop if dealing with less spectra or a list of unique spectra wavelength ranges
  for i = 0,54 do begin
    ;get the object spectrum; called 'object' here
    fits_read, spectrafiles[i], object, header
    wavelengthRange = (sxpar(header,'CRVAL1') - (sxpar(header,'CRPIX1') - 1d) * sxpar(header,'CDELT1')) + [0d,sxpar(header,'CDELT1')*(sxpar(header,'NAXIS1') - 1d)]
     ;converting from log to linear if needed
     ;wavelengthRange = 10^(wavelengthRange)
    log_rebin, wavelengthRange, object, logObject, logWavelengths, VELSCALE=velocityScale


    ;get the error array; called 'error' here
    fits_read, noisefiles[i], error, header
    ;Convert the ".ivar" file values into the raw noise values ("ivar" is the inverse variance)
    error = sqrt(abs(double(error)))
    error = double(1./error)
    log_rebin, wavelengthRange, error, logError, logWavelengths


    ;get the wavelength range of the template library
    fits_read, templatefiles[0], template, header
      wavelengthRange = (sxpar(header,'CRVAL1') - (sxpar(header,'CRPIX1') - 1d) * sxpar(header,'CDELT1')) + [0d,sxpar(header,'CDELT1')*(sxpar(header,'NAXIS1') - 1d)]
      log_rebin, wavelengthRange, template, logTemplate, logTemplateWavelengths, VELSCALE=velocityScale
    ;get the size of the template spectra library
    templates = dblarr(n_elements(logTemplate), ntemplates)
    ;make a 2D array that is the number of template files by the number of pixels (called 'template' 5 lines above) in each file's spectrum
    for j = 0,ntemplates-1 do begin
      fits_read, templatefiles[j], template, header
      wavelengthRange = (sxpar(header,'CRVAL1') - (sxpar(header,'CRPIX1') - 1d) * sxpar(header,'CDELT1')) + [0d,sxpar(header,'CDELT1')*(sxpar(header,'NAXIS1') - 1d)]
      log_rebin, wavelengthRange, template, logTemplate, logTemplateWavelengths, VELSCALE=velocityScale
      templates[*,j] = logTemplate
    endfor


    output = dblarr(2,n_elements(object))
 
    logWavelengthMask = alog(wavelengthMask)
    mask = dblarr(n_elements(logObject))
    listSize = size(logWavelengthMask)

    for i = 0, listSize[2] - 1 do begin
      inrange = (logWavelengths gt logWavelengthMask[0,i]) and (logWavelengths lt logWavelengthMask[1,i])
      mask = inrange or mask
    endfor

    pixelmask = where(not mask)
    ;pixelmask = indgen(n_elements(logObject))
    goodPixels = pixelmask

    velocityOffset = (logTemplateWavelengths[0] - logWavelengths[0]) * 299792.458d; speed of light, km/s



    ;Changeable parameters, defaults are given on the right.
    multiDegree = 7;7
    moments = 4;2
    degs = 1;-1
 
 
 
    
    ;Example loop for iterating over velocities if none are known for a particular spectrum
    for i = 0,799, 100 do begin; Iterate over velocities to find a suitable input
      input_vel = float(50+i) ;start_vel[s] ;km/s
      input_disp = float(30);start_sig[s]; km/s
    ;start = [veloc[s], 30.] ;starting guess values for velocity and velocity dispersion from list of previously measure values
    start = [input_vel, input_disp] ;starting guess from for loop (for first time use)

    ppxf, templates, logObject, logError, velocityScale, start, solution, /CLEAN, BESTFIT=fit, $
      GOODPIXELS=goodPixels, VSYST=velocityOffset, BIAS=0, WEIGHTS=weights, $
      MOMENTS=moments, MDEGREE=multiDegree, DEGREE=degs, /QUIET, /PLOT

    ;print, 'Formal errors:    dV    dsigma       dh3       dh4'
    ;print, error*sqrt(solution[6]), FORMAT='(10x,2f10.1,2f10.3)'

    ;Optional resampling routine. This adds artificial noise to the object spectra and re-fits them with pPXF.
    for i = 0, resamplings - 1 do begin
      noise = logError * randomn(seed, n_elements(logError))
      ;Add artificial noise to resampling routine
      sample = fit + noise
      goodPixels = pixelmask

      print,'Resample '+strtrim(i+1,2)+''
      ppxf, templates, sample, usederror, velocityScale, start, solution, $
        /CLEAN, BESTFIT=sampleFit, GOODPIXELS=goodPixels, MOMENTS=moments, $
        MDEGREE=multiDegree, VSYST=velocityOffset, BIAS=0, WEIGHTS=weights, DEGREE=degs, /QUIET, $
        /PLOT

    endfor

    ;print,''
    ;print, 'Velocity (km/s): ', solution[0]
    ;print, 'Dispersion (km/s): ', solution[1]
    ;print, 'Chi2/DOF: ', solution[6]
    ;print,''

    ;endfor

    printf,1,''

  endfor

  CLOSE, /ALL



end
