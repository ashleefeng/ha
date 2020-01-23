
pro smb_peakmovie_maker_3alex, run, text_ID, color_number_input

	;3color blue-green-red 3ALEX
	color_number = 3
	; Custumizing parameters
	movie_width = 11		;512 image-> 15, 256 image-> 11

	loadct, 5
	device, decomposed=0

	COMMON colors, R_ORIG, G_ORIG, B_ORIG, R_CURR, G_CURR, B_CURR

	if N_PARAMS() eq 0 then begin		;run arrgument is empty.
		run = DIALOG_PICKFILE(PATH='c:\user\tir', TITLE='Select pma File', /READ, FILTER = '*.pma')
		xdisplayFile, '', TEXT=(run + "is selected to make a movie."), RETURN_ID=display_ID, WTEXT=text_ID
		run = strmid(run, 0, strlen(run) - 4)
	endif
	
	if N_PARAMS() eq 1 then begin
		xdisplayFile, '', TEXT=(run + ".pma is selected to make a movie."), RETURN_ID=display_ID, WTEXT=text_ID
	endif
	
	if N_PARAMS() eq 2 then begin
		WIDGET_CONTROL, text_ID, SET_VALUE=(run + ".pma is selected to make a movie."), /APPEND, /SHOW
	endif
	
	if N_PARAMS() eq 3 then begin
		WIDGET_CONTROL, text_ID, SET_VALUE=(run + ".pma is selected to make a movie."), /APPEND, /SHOW
		color_number = color_number_input
	endif
	

	; figure out size + allocate appropriately
	close, 1				; make sure unit 1 is closed
	openr, 1, run + ".pma"

	file_infomation = FSTAT(1)

	film_width = fix(1)
	film_height = fix(1)
	readu, 1, film_width
	readu, 1, film_height

	film_time_length = long(       long64(file_infomation.SIZE -long64(4))/(long64(film_width)*long64(film_height))     )

	WIDGET_CONTROL, text_ID, SET_VALUE=("film width, height, time_length : " + STRING(film_width) + STRING(film_height) + STRING(film_time_length)), /APPEND, /SHOW


	; load the locations of the peaks
	;answer=gui_prompt('How many color channel on the display (2 or 3):', title='2 or 3 color')
	;colornumber = total(long(answer))
	
	close, 2
	openr, 2, run + ".3color_3alex_pks"
	
	GoodLocations_x = intarr(4000)
	GoodLocations_y = intarr(4000)
	Background_blue = dblarr(4000)
	Background_green = dblarr(4000)
	Background_red = dblarr(4000)
	NumberofGoodLocations = fix(1)
	xx = fix(1)
	yy = fix(1)
	back_blue = double(1)
	back_green = double(1)
	back_red = double(1)
	color_change_check='n'
	readf, 2, NumberofGoodLocations
	for i = 0, NumberofGoodLocations - 1 do begin
    	readf, 2, indexdummy, xx, yy, back_blue, back_green, back_red
    	GoodLocations_x(i) = xx
    	GoodLocations_y(i) = yy
    	Background_blue(i) = back_blue
    	Background_green(i) = back_green
    	Background_red(i) = back_red
	endfor
	readf, 2, color_change_check
	close, 2

  WIDGET_CONTROL, text_ID, SET_VALUE=("Color change? : " + color_change_check), /APPEND, /SHOW

	WIDGET_CONTROL, text_ID, SET_VALUE=(STRING(NumberofGoodLocations/color_number) + " peaks were found in file " + run + ".pma"), /APPEND,  /SHOW

	close, 2

	close, 3				; make sure unit 3 is closed
	openw, 3, run + ".3color_3alex_movies"

	;answer=gui_prompt('movie width (recommanded=15 pixel) :', title='Movie width')
	;movie_width=round(total(long(answer)))

	half_movie_width = (movie_width-1)/2
	movie_width_total =movie_width*NumberofGoodLocations
	writeu, 3, fix(movie_width_total)
	writeu, 3, fix(movie_width)

	frame = bytarr(film_width, film_height)
	temp= fltarr(film_width, film_height)

	;if color_change_check eq 'y' then begin
	;	film_time_length = film_time_length-1
	;	readu, 1, frame
	;endif

	peaks = bytarr(movie_width_total, movie_width)

	WIDGET_CONTROL, text_ID, SET_VALUE="Start.", /APPEND,  /SHOW
	for t= 0,film_time_length-1 do begin
		if (t mod 100) eq 0 then begin
			WIDGET_CONTROL, text_ID, SET_VALUE=("Movie Working on : " + STRING(t) +  "/" + STRING(film_time_length)), /APPEND, /SHOW
		endif
		readu, 1, frame
		for j = 0, NumberofGoodLocations - 1 do begin
		  if (t mod 3) eq 0  then begin
        temp = float(frame((GoodLocations_x(j)-half_movie_width):(GoodLocations_x(j)+half_movie_width), (GoodLocations_y(j)-half_movie_width):(GoodLocations_y(j)+half_movie_width))) - Background_blue(j)
      endif
			if (t mod 3) eq 1  then begin
				temp = float(frame((GoodLocations_x(j)-half_movie_width):(GoodLocations_x(j)+half_movie_width), (GoodLocations_y(j)-half_movie_width):(GoodLocations_y(j)+half_movie_width))) - Background_green(j)
			endif
			if (t mod 3) eq 2  then begin
				temp = float(frame((GoodLocations_x(j)-half_movie_width):(GoodLocations_x(j)+half_movie_width), (GoodLocations_y(j)-half_movie_width):(GoodLocations_y(j)+half_movie_width))) - Background_red(j)
			endif
			peaks((j*movie_width):((j+1)*movie_width-1), *) = byte(temp>0)
		endfor
		writeu, 3, peaks
	endfor

	close, 1
	close, 3

	WIDGET_CONTROL, text_ID, SET_VALUE="Movie maker for " + run +".pma file", /APPEND,  /SHOW
	WIDGET_CONTROL, text_ID, SET_VALUE="Done. ", /APPEND,  /SHOW

end

