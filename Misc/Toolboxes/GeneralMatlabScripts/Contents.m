% ASTRONOMICAL - General utilities
%                                 By : Eran O. Ofek
%                                 Version: May 2006
%
% List of MATLAB programs in the general package
%
% and_nan   - Logical function "and" for NaNs Similar to "and" logical
%             function, but NaNs are regarded as no information.
% bc_a      - Bias correction and acceleration for
%             Bootstrap and Jackknife estimetes of the
%             confidence interval.
%             The bc_a bias coorection and acceleration can be used
%             to estimate the bias corrected and accelerated
%             confiedence interval.
% bin_sear  - Binary search of nearest value in a given sorted vector.
% binsear_f - Binary search on monotonic function
%             (Given monotonic function y=f(x) and y, search for x).
% bivar_gauss - Return the value of a normalized bivariate gaussian
%             in a list of points.
% bootstrap_std - Given an estimator (given by a function)
%             calculate the Bootstrap StD (and bias) for this estimator.
% calc_limitmag - Calculate eproximate limiting magnitude.
%             Given telescope parameters, and S/N calculate the
%             limiting magnitude for a point source.
% cel_coo_rnd- Generate random coordinates on the celestial sphere.
%             The program can applay matrix rotation reject
%             coordinates, and generate non uniform coordinates.
% cell_stat - Given a list of x,y coordinates (with optional
%             property columns), count the number of points in
%             each cell in the xy space and some statistical
%             properties (the mean, median,... ignore nans).
% convert_energy -Convert between different energy units.
% convert_units - Unit conversion function
%             Given an input and output strings
%             containing units - return the conversion
%             multiplication factor needed for converting
%             the input units to the output units.
% delete_ind- Delete a column/s or row/s from a specific position in a matrix.
% differr   - Calculate differences and errors.
%             If two input arguments V1, V2 are given then
%             the program returns:
%             [V1(:,1)-V2(:,1), sqrt(V1(:,2).^2+V2(:,2).^2)]
%             If one input argument is given the the program
%             differentiate successive values.
% dist_p2line-Calculate the minimum distance in a 2-d space
%             between a line and a point.
% ensemble_from_error- Generate normal distributed random
%             ensamble realization given a mean value
%             and the one (or two) sided standard
%             deviation.
% equipot   - Calculate two body equipotanials map.
% err_cl    - Calculate right/left confidence boundary to a
%             given numerical distribution.
% find_str_cell - Search for a (sub)string within cell array
%             of strings and return the indicies of the
%             cell-array elements that contain the (sub)string.
% fit_circle - Fit points, on a plane or a sphere, to
%             a circle. Calculate the best fit radius and
%             the center of the circle.
% fit2d     - Find 2D geometrical transformation
%             between two sets of x/y coordinates.
%             The first set is the reference set
%             and the second is the data.
%             Fit of the form:
%             Y' = dy + ay0(y-y0)   + by0(x-x0)   +
%                       ay1(y-y0)^2 + by1(x-x0)^2 +
%                       cy(x-x0)(y-y0)
%             X' = dx + ax0(x-x0)   + bx0(y-y0)   +
%                       ax1(x-x0)^2 + bx1(y-y0)^2 +
%                       cx(x-x0)(y-y0)
%             Where x/y are the reference coordinates
%             and X'/Y' are the data coordinates.
%             x0/y0 are coordinates defined by the user.
%             a0/a1/b0/b1/c/d are the model parameters.
% for_each_file - Given a file name containing list of files,
%             load each file into a matrix and execute
%             a function with the loaded matrix as a paramter.
% fpf       - Easy to use fprintf.
%             Similar to fprintf, but open and close the file
%             automaticaly. In case that format string is
%             not given then the function select a format
%             automaticaly based on number precision.
% fresnelc  - Fresnel cosine function, return: cos(0.5*pi*T^2).
% fresnels  - Fresnel sine function, return: sin(0.5*pi*T^2).
% get_constant - Get astronomical/physical constant from list of constants.
% html_table- Given a matlab matrix or cell array
%             create an html page with an html table.
% image2avi - Create avi file from list of images
% integral_percentile - Given a function (X,Y), calculate the limits of
%             the integal which contains a given percentile of the
%             total integral.
% insert_ind- Insert a column/s or row/s to a specific position in a matrix.
% jackknife - Given an estimator (given by a function)
%             calculate the Jackknife StD and the first order
%             Quenouille-Tukey jacknife bias for
%             this estimator.
% latex_table - Create latex table from a data given in the cell array.
% linemindist_t - Given two 2D linear trajectories
%             calculate the time in which the distance
%             between the two trajectories is minimal
%             and the distance at that time.
%             Each line is defined by:
%             x_i = A_i + B_i*(t - T_i)
%             y_i = D_i + E_i*(t - T_t)
% load2cell - Read table into cell array.
%             This is useful whenever there is text in the table.  
%             Numbers are stored as numbers and empty 'cells' are
%             filled with NaNs.
% mat2vec   - convert matrix to vector.
% match_lists_index- Match lists according to index.
% match_lists_position- Match lists according to 2-D position.
% match_lists_time -Match lists according to 1-D position (e.g., time).
% max_likelihood- Given a numerical distribution and list
%             of "events". Calculates the maximum likelihood of events.
%             Calculate the ML ratio test.
%             Using Monte-Carlo simulation of the parent
%             numerical distribution, calculate the ML
%             probability distribution.
% max2d     - 2d maximum function. Return the maximum value
%             and index in a 2d matrix.
% maxnd     - Return the global maximum of a N-D matrix and its indices.
% meannd    - Return the global mean of a N-D matrix.
% mediannd  - Return the global median of a N-D matrix.
% merge_by_coo- merge two catalogs by coordinates.
%             assume coordinates are in radians.
%             reject entries that do not have a match.
% mgrep     - grep-like utility for MATLAB.
%             Search for substrings within text file.
% min2d     - 2d minimum function. Return the minimum value
%             and index in a 2d matrix.
% minnd     - Return the global minimum of a N-D matrix and its indices.
% min_sear  - searching for value position in a given vector by
%             minimazation.
% mode      - Estimate the mode value in a given vector
%             (The most probable value).
% or_nan    - Logical function "or" for NaNs Similar to "or" logical
%             function, but NaNs are regarded as no information.
% poissconf - Given the number of observed events, assuming
%             Poisson statistics, calculate the two sided upper
%             and lower confidence intervals.
% prob2find_inr - Given an object density (number per unit area),
%             and a distance from a point, calculate the
%             the probability to find an object within radius R
%             from the point (assuming Poisson statistics).
% prob_ci   - Calculate 1,2,3\sigma (and others) confidence
%             interval from a given numerical probability
%             function.
% rand_circle- Generate random number equally distributed
%             inside a unit circle.
% randinpolygon - Generate random positions inside (or on boundry)
%             a polygon, on a plane or a sphere.
% randgen   - General random numbers generator.
%             Generate vector of random numbers from a given
%             numerical distribution.
% rand_ps   - Generate a random realization of a time series
%             with a given power spectra (e.g., power-law) and
%             optional gaussian measuments errors.
% rangend   - Return the global Range of a N-D matrix.
% realhist  - calculate histogram on a pregiven range.
% remove_cell_element - Remove a list of indices from a cell vector.
% rotm      - return the rotation matrix about the X/Y/Z axis.
% stdnd     - Return the global StD of a N-D matrix.
% sn_calc   - Calculate eproximate Signal to Noise (S/N) ratio.
%             Given telescope parameters, calculate the S/N of
%             a point source.
% sn_diag   - Sum vs. Number diagram. Total population up to
%             number vs. number.
% stellar_imf - Return the stellar initial mass function in a
%             given mass range.
% symerror  - Given a symbolic expression and the
%             variables of the expression, calculate
%             the symbolic error function of the
%             expression, with respect to the variables.
%             The new expression contains an error
%             variable named D_"original_var" for each
%             variable respectively.
%             Use: char(vectorize(ErrExp)) to convert
%             the symbolic expression to a vectorized
%             string that can be evaluated using eval.
%             The function also returns a cell array of
%             the new error variables
%             (e.g.,  D_"original_var").
% wc        - Call Unix wc (word count) command.
% wmean     - wighted mean.
% xyxymatch - Match two lists by 2-d coordinates (spherical
%             or planner). The program indicate about possible
%             confusion in the matching (i.e., one to many...).
%
