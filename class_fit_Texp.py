
####  This class is used to fit a polynomial to the light curve to look for the explosion time (time of zero flux) 
## Make sure to convert magnitude measurements to  flux measurements


from iminuit.util import propagate


from library_funciton_fitTexp import *

class Fit_t_exp_flux( object ):

    #TODO: fix the chi2on the full data

    '''
        This class is defined to fit the early LC of Infant SNe II in order to estimate the explosion parameters
        The object given to the class shoul be "time (from detection), mag, e_mag"

    '''

    def __init__(self, table, min_fit, max_fit, mini = -10, maxi = 7, zero_point=None):
        '''
        Parameters
        ----------
        table contraining the time from first marshal detection, the flux and the error on the flux 

        min_fit    [float] minimal bound of the fit
        max_fit    [float] maximal bound of the fit
        mini       [flaot] minimal bound of data considered
        maxi       [flaot] maximal bound of data considered for the final chi_2 fit 
        zero_point [float] MAGNITUDE of what you will consider as your "zero" for fitting purposes. 
                        This can be the levels of non detections or of the progenitor




        '''
        self._set_max_fit_val(max_fit)
        self._time_zone_selected = False
        # self._filter_name   = {'ZTF_g':'darkolivegreen','ZTF_r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        # self.set_filter(filter)

        # full data to comapre the chi_2 to
        self.set_time_zone_interest(mini, maxi)
        self.set_time_from_asfd()
        self.set_flux(zero_point)
        self.set_e_flux()


        # reduced data to perform the fit onto
        self.set_time_zone_fit(min_fit, max_fit)
        # if len(self.fittable) <

        self.set_time_from_asfd_fit()
        self.set_flux_fit(zero_point)
        self.set_e_flux_fit()


        self._minifitplot = mini
        self._maxifitplot = maxi

        

    #-----------#
    #  SETTERS  #
    #-----------#

    def _set_max_fit_val(self, max_fit):
        '''
        '''
        self._max_fit = max_fit

    def set_filter(self,filter):
        '''
        '''
        self.filter = filter

    def set_phot_table(self, table):
        '''

        '''
        self.table = table

    def set_time_zone_interest(self, mini = -2, maxi = 0.5):
        '''
        the early light curve we're interested in:
        taking from 10 days prior to FD until a week from FD
        '''
        self.table = self.table[(self.table['t_fromfd']>=mini)&(self.table['t_fromfd']<=maxi)]

    def set_time_from_asfd(self):
        '''
        '''
        self.t = self.table['t_fromfd']

    def set_flux(self, zero_point):
        '''
        '''


        self.f = to_flux(self.table['mag'])

        if zero_point is not None:
            self._zero_flux = to_flux(zero_point)

            self.f = self.f - self._zero_flux

    def set_e_flux(self):
        '''
        '''

        self.e_f = error_onflux_from_mag(self.table['mag'],self.table['magerr'])




    def set_time_zone_fit(self, min_fit , max_fit ):
        '''

        '''

        self.fittable = self.table[(self.table['t_fromfd']>=min_fit)&(self.table['t_fromfd']<=max_fit)]

        self._time_zone_selected = True

    def set_time_from_asfd_fit(self):
        '''
        '''
        self.t_fit = self.fittable['t_fromfd']

    def set_flux_fit(self,zero_point):
        '''
        '''

        self.f_fit = to_flux(self.fittable['mag'])
        if zero_point is not None:
            self._zero_flux = to_flux(zero_point)
            
            self.f = self.f - self._zero_flux

    def set_e_flux_fit(self):
        '''
        '''
        self.e_f_fit = error_onflux_from_mag(self.fittable['mag'],self.fittable['magerr'])


    ############
    # GETTERS  #
    ############


    # def get_fit_roi_length(self):
    #     '''
    #     This function gets the length of the region of ionterest of fit (early rise)
    #     '''
    #     if self._time_zone_selected == True:
    #         len_roi = len(self.fittable[self.fittable['t_fromfd']>=0])
    #         return len_roi
    #     else:
    #         print('You need to select a region of interest got the fit')

    
    def get_chi_2(self, param):
        '''
        This function returns the Chi squared of the function fitted to the data. 

        parameters
        ----------
        self

        returns
        -------

        ''' 
        #a , T_exp , n  = param 
        model_ = rise_typeII(self.t_fit,*param)
        pull_  = (( self.f_fit - model_ )/ self.e_f_fit )**2
        chi_sq = np.sum(pull_)
        dof    = len(self.f_fit)-len(param)
        if dof == 0: #avoid division by zero? 
            dof = 1

        return  chi_sq / dof
    
    def _get_chi2minuit(self, a, t_exp, n):
        '''

        '''

        param = [a,t_exp,n]

        return self.get_chi_2(param)



    def get_explosion_t(self):
        '''
        '''
        self.t_exp  = self.minuit_output.params[1].value
        self.dt_exp = self.minuit_output.params[1].error
        return self.t_exp, self.dt_exp
        

    def get_fit_index(self):
        '''
        '''
        self.n      = self.minuit_output.params[2].value
        self.dt_n   = self.minuit_output.params[2].error
        return self.n, self.dt_n




    #######################
    #        FITTER       #
    #######################
    
    def set_minuit(self, guess, boundaries=None, fixed=None, step_sizes = None, print_level=0):
        '''
        
        '''

        a, t_exp, n     = guess
        self._paramname = "a,t_exp,n".split(",")
        
        local_ = locals() # dict with the local variables


        self.minuit_kwargs = {}
        for p in self._paramname:
            self.minuit_kwargs[p] = local_[p]


        # initialise the optimiser: we are using our own cost function defined as a reduced chi_2 and we feed it the guesses
        self.minuit = Minuit(self._get_chi2minuit, **self.minuit_kwargs)

        # since we are using a least square type of cost funtion, we need to define the error level definition here

        self.minuit.errordef = Minuit.LEAST_SQUARES

        self.minuit.print_level = print_level

        # we then set all the parameters supposedly fixed, their limits and step sizes 
        if boundaries is not None:
            for _ in boundaries.keys():
                self.minuit.limits[_] = boundaries[_]

        if fixed is not None:
            for _ in fixed:
                self.minuit.fixed[_] = True

        if step_sizes is not None: 
            for _ in step_sizes.keys():
                self.minuit.errors[_] = step_sizes[_]





        # self.minuit_kwargs = {}

        # for p in self._paramname:
        #     self.minuit_kwargs[p] = local_[p]
            
        # if boundaries is not None:
        #     for k,v in boundaries.items():
        #         self.minuit_kwargs["limit_"+k] = v
                
        # if fixed is not None:
        #     for k in fixed:
        #         self.minuit_kwargs["fix_"+k] = True
        
        # self.minuit = Minuit(self._get_chi2minuit, errordef=errordef, error_a = 5e-9 , error_t_exp = 0.1, error_n = 0.05,
        #                             print_level=print_level, **self.minuit_kwargs)
        

    def fit_minuit(self, guess, boundaries=None, fixed=None, step_sizes=None):
        
        """ 
    
        """
        
        self.set_minuit(guess=guess, boundaries=boundaries, fixed=fixed, step_sizes = step_sizes )
        self.minuit_output = self.minuit.migrad()
        print(self.minuit_output)

    def error_propag(self, min_interp = None, max_interp = None, nump = 500 , returnop = False ) :
        '''
        this function propagates the error on the minuit fit and returns the estimated errors 

        parameters
        ----------
        min_interp [float]
        max_interp [float]
        nump       [float] number of points to interpolate, pd 500
        
        '''
        if (min_interp, max_interp) == (None, None):
            min_interp = min(self.t_fit)
            max_interp = max(self.t_fit)

        self.t_res = np.linspace(min_interp, max_interp, nump)

        self.fit_mag_res , self.fit_mag_cov_res  = propagate(lambda param: rise_typeII(self.t_res, *param), self.minuit.values, self.minuit.covariance)
        
        self.fit_e_mag_res = np.diag(self.fit_mag_cov_res)**0.5

        if returnop is True:
            return self.t_res, self.fit_mag_res, self.fit_e_mag_res  
    
    ##################
    #   PLOTTING     #
    ##################



    def plot_fit(self, add_plot = None, scale = None):
        '''
        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_ = np.linspace(self._minifitplot, self._maxifitplot, 1000)
        
        _params = [self.minuit_output.params[x].value for x in range(3)]
        _errors = [self.minuit_output.params[x].error for x in range(3)]



        # t_ , fit, e_fit = self.error_propag(returnop=True)
        
        tx_string = f't_exp = {_params[1]:.3f}±{_errors[1]:.3f} d '
        n_string  = f'n = {_params[2]:.3f}±{_errors[2]:.3f} '

        # chi2_      = self.get_redchi_2_gen(_params)
        # chi_string = f'Chi_2 early LC = {chi2_:.3f}'

        # chi2_fit       = self.get_redchi_2_fit(_params)
        # chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        # fit_interval_string = f'd = {self._max_fit}' 




        if add_plot == None:
            plt.figure()
            plt.errorbar(self.t, self.f, self.e_f, fmt = 'o', ms=3.5 , color = 'grey' )
            plt.plot(x_ , rise_typeII(x_, *_params), ls = '--' , label= n_string , color = _fitcolor )

            # plt.plot(x_[0] , rise_typeII(x_, *_params)[0], color = 'white',label= chi_string )
            # plt.plot(x_[1] , rise_typeII(x_,*_params)[1], color = 'white',label= chi_fit_string )


            plt.axvline(_params[1], color = _fitcolor ,ls=':',alpha = 0.5, label = tx_string )
            plt.axvspan(_params[1] - _errors[1], _params[1] + _errors[1], color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from ASFD [days]', size = 15)
            plt.ylabel('Flux - Fo [Jy]', size = 15)

            if scale is not None:
                plt.yscale("log")

            plt.legend()
        
        # else: 
        #     add_plot.plot(0,0, ls = '--' , label= fit_interval_string , color = 'white' )
        #     add_plot.errorbar(self.t, self.f, self.e_f, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
        #     add_plot.plot(x_ , rise_typeII(x_, a_, t_exp_, n_), ls = '--' , label= n_string , color = _fitcolor )

        #     add_plot.plot(x_[0] , rise_typeII(x_, a_, t_exp_, n_)[0], color = 'white',label= chi_string )
        #     add_plot.plot(x_[1] , rise_typeII(x_, a_, t_exp_, n_)[1], color = 'white',label= chi_fit_string )


        #     add_plot.axvline(t_exp_, color = _fitcolor ,ls=':',alpha = 0.5, label = tx_string )
        #     # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


        #     plt.xlabel('Time from ASFD [days]', size = 15)
        #     plt.ylabel('Flux [Jy]', size = 15)

        #     plt.legend(fontsize = 10)

        
        
        

        
        
 


    # def get_redchi_2_gen(self, param):
    #     '''
    
    #     parameters
    #     ----------
    #     self

    #     returns
    #     -------

    #     ''' 
    
    #     #a , T_exp , n  = param 
    #     model_ = rise_typeII(self.t,*param)
    #     pull_  = (( self.f - model_ )/ self.e_f )**2
    #     chi_sq = np.sum(pull_)
    #     dof    = len(self.f) - len(param)

    #     return  chi_sq/dof 


    # def get_redchi_2_fit(self, param):
    #     '''
    #     parameters
    #     ----------
    #     self

    #     returns
    #     -------
    #     ''' 
    #     #a , T_exp , n  = param 
    #     model_ = rise_typeII(self.t_fit,*param)
    #     pull_  = (( self.f_fit - model_ )/ self.e_f_fit )**2
    #     chi_sq = np.sum(pull_)
    #     dof    = len(self.f_fit)-len(param)

    #     return  chi_sq/dof




















        
'''
##############################################################################################################






THIS IS THE FITTER IN MAGNITUDE SPACE 




##############################################################################################################

'''


class Fit_t_exp_mag( object ):

    #TODO: fix the chi2on the full data

    '''
        This class is defined to fit the early LC of Infant SNe II in order to estimate the explosion parameters
        The object given to the class shoul be "time (from detection), mag, e_mag"

    '''

    def __init__(self, table, min_fit, max_fit, mini = -10, maxi = 7):
        '''
        Parameters
        ----------
        table contraining the time from first marshal detection, the flux and the error on the flux 

        min_fit [float] minimal bound of the fit
        max_fit [float] maximal bound of the fit
        mini    [flaot] minimal bound of data considered
        maxi    [flaot] maximal bound of data considered for the final chi_2 fit 




        '''
        self._set_max_fit_val(max_fit)
        self._time_zone_selected = False
        # self._filter_name   = {'ZTF_g':'darkolivegreen','ZTF_r':'orangered', 'ZTF_i':'gold'}

        self.set_phot_table(table)
        # self.set_filter(filter)

        # full data to comapre the chi_2 to
        self.set_time_zone_interest(mini, maxi)
        self.set_time_from_asfd()
        self.set_mag()
        self.set_e_mag()


        # reduced data to perform the fit onto
        self.set_time_zone_fit(min_fit, max_fit)
        # if len(self.fittable) <

        self.set_time_from_asfd_fit()
        self.set_mag_fit()
        self.set_e_mag_fit()


        self._minifitplot = mini
        self._maxifitplot = maxi

        

    #-----------#
    #  SETTERS  #
    #-----------#

    def _set_max_fit_val(self, max_fit):
        '''
        '''
        self._max_fit = max_fit

    def set_phot_table(self, table):
        '''

        '''
        self.table = table

    def set_time_zone_interest(self, mini = -3, maxi = 0.5):
        '''
        the early light curve we're interested in:
        taking from 10 days prior to FD until a week from FD
        '''
        self.table = self.table[(self.table['t_fromfd']>=mini)&(self.table['t_fromfd']<=maxi)]

    def set_time_from_asfd(self):
        '''
        '''
        self.t = self.table['t_fromfd']

    def set_mag(self):
        '''
        '''


        self.m = self.table['mag']

    def set_e_mag(self):
        '''
        '''

        self.e_m = self.table['magerr']




    def set_time_zone_fit(self, min_fit , max_fit ):
        '''

        '''

        self.fittable = self.table[(self.table['t_fromfd']>=min_fit)&(self.table['t_fromfd']<=max_fit)]

        self._time_zone_selected = True

    def set_time_from_asfd_fit(self):
        '''
        '''
        self.t_fit = self.fittable['t_fromfd']

    def set_mag_fit(self):
        '''
        '''

        self.m_fit = self.fittable['mag']

    def set_e_mag_fit(self):
        '''
        '''
        self.e_m_fit = self.fittable['magerr']


    ############
    # GETTERS  #
    ############


    # def get_fit_roi_length(self):
    #     '''
    #     This function gets the length of the region of ionterest of fit (early rise)
    #     '''
    #     if self._time_zone_selected == True:
    #         len_roi = len(self.fittable[self.fittable['t_fromfd']>=0])
    #         return len_roi
    #     else:
    #         print('You need to select a region of interest got the fit')

    
    def get_chi_2(self, param):
        '''
        This function returns the Chi squared of the function fitted to the data. 

        parameters
        ----------
        self

        returns
        -------

        ''' 
        #a , T_exp , n  = param 
        model_ = rise_typeII_mag(self.t_fit,*param)
        pull_  = (( self.m_fit - model_ )/ self.e_m_fit )**2
        chi_sq = np.sum(pull_)
        dof    = len(self.m_fit)-len(param)
        if dof == 0: #avoid division by zero? 
            dof = 1

        return  chi_sq / dof
    
    def _get_chi2minuit(self, a, t_exp, n):
        '''

        '''

        param = [a,t_exp,n]

        return self.get_chi_2(param)



    def get_explosion_t(self):
        '''
        '''
        self.t_exp  = self.minuit_output.params[1].value
        self.dt_exp = self.minuit_output.params[1].error
        return self.t_exp, self.dt_exp
        

    def get_fit_index(self):
        '''
        '''
        self.n      = self.minuit_output.params[2].value
        self.dt_n   = self.minuit_output.params[2].error
        return self.n, self.dt_n




    #######################
    #        FITTER       #
    #######################
    
    def set_minuit(self, guess, boundaries=None, fixed=None, step_sizes = None, print_level=0):
        '''
        
        '''

        a, t_exp, n     = guess
        self._paramname = "a,t_exp,n".split(",")
        
        local_ = locals() # dict with the local variables


        self.minuit_kwargs = {}
        for p in self._paramname:
            self.minuit_kwargs[p] = local_[p]


        # initialise the optimiser: we are using our own cost function defined as a reduced chi_2 and we feed it the guesses
        self.minuit = Minuit(self._get_chi2minuit, **self.minuit_kwargs)

        # since we are using a least square type of cost funtion, we need to define the error level definition here

        self.minuit.errordef = Minuit.LEAST_SQUARES

        self.minuit.print_level = print_level

        # we then set all the parameters supposedly fixed, their limits and step sizes 
        if boundaries is not None:
            for _ in boundaries.keys():
                self.minuit.limits[_] = boundaries[_]

        if fixed is not None:
            for _ in fixed:
                self.minuit.fixed[_] = True

        if step_sizes is not None: 
            for _ in step_sizes.keys():
                self.minuit.errors[_] = step_sizes[_]





        # self.minuit_kwargs = {}

        # for p in self._paramname:
        #     self.minuit_kwargs[p] = local_[p]
            
        # if boundaries is not None:
        #     for k,v in boundaries.items():
        #         self.minuit_kwargs["limit_"+k] = v
                
        # if fixed is not None:
        #     for k in fixed:
        #         self.minuit_kwargs["fix_"+k] = True
        
        # self.minuit = Minuit(self._get_chi2minuit, errordef=errordef, error_a = 5e-9 , error_t_exp = 0.1, error_n = 0.05,
        #                             print_level=print_level, **self.minuit_kwargs)
        

    def fit_minuit(self, guess, boundaries=None, fixed=None, step_sizes=None):
        
        """ 
    
        """
        
        self.set_minuit(guess=guess, boundaries=boundaries, fixed=fixed, step_sizes = step_sizes )
        self.minuit_output = self.minuit.migrad()
        print(self.minuit_output)

    def error_propag(self, min_interp = None, max_interp = None, nump = 500 , returnop = False ) :
        '''
        this function propagates the error on the minuit fit and returns the estimated errors 

        parameters
        ----------
        min_interp [float]
        max_interp [float]
        nump       [float] number of points to interpolate, pd 500
        
        '''
        if (min_interp, max_interp) == (None, None):
            min_interp = min(self.t_fit)
            max_interp = max(self.t_fit)

        self.t_res = np.linspace(min_interp, max_interp, nump)

        self.fit_mag_res , self.fit_mag_cov_res  = propagate(lambda param: rise_typeII_mag(self.t_res, *param), self.minuit.values, self.minuit.covariance)
        
        self.fit_e_mag_res = np.diag(self.fit_mag_cov_res)**0.5

        if returnop is True:
            return self.t_res, self.fit_mag_res, self.fit_e_mag_res  
    
    ##################
    #   PLOTTING     #
    ##################



    def plot_fit(self, add_plot = None):
        '''
        '''

        _colors  = ['blue', 'teal', 'lightblue', 'dimgrey','rebeccapurple']
        _fitcolor = random.choice(_colors)


        x_ = np.linspace(self._minifitplot, self._maxifitplot, 1000)
        
        _params = [self.minuit_output.params[x].value for x in range(3)]
        _errors = [self.minuit_output.params[x].error for x in range(3)]



        # t_ , fit, e_fit = self.error_propag(returnop=True)
        
        tx_string = f't_exp = {_params[1]:.3f}±{_errors[1]:.3f} d '
        n_string  = f'n = {_params[2]:.3f}±{_errors[2]:.3f} '

        # chi2_      = self.get_redchi_2_gen(_params)
        # chi_string = f'Chi_2 early LC = {chi2_:.3f}'

        # chi2_fit       = self.get_redchi_2_fit(_params)
        # chi_fit_string = f'Chi_2 fit = {chi2_fit:.3f}'

        # fit_interval_string = f'd = {self._max_fit}' 




        if add_plot == None:
            plt.figure()
            plt.errorbar(self.t, self.m, self.e_m, fmt = 'o', ms=3.5 , color = 'grey' )
            plt.plot(x_ , rise_typeII_mag(x_, *_params), ls = '--' , label= n_string , color = _fitcolor )

            # plt.plot(x_[0] , rise_typeII(x_, *_params)[0], color = 'white',label= chi_string )
            # plt.plot(x_[1] , rise_typeII(x_,*_params)[1], color = 'white',label= chi_fit_string )


            plt.axvline(_params[1], color = _fitcolor ,ls=':',alpha = 0.5, label = tx_string )
            # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


            plt.xlabel('Time from ASFD [days]', size = 15)
            plt.ylabel('Flux [Jy]', size = 15)
            plt.gca().invert_yaxis()
            plt.legend()
        
        # else: 
        #     add_plot.plot(0,0, ls = '--' , label= fit_interval_string , color = 'white' )
        #     add_plot.errorbar(self.t, self.f, self.e_f, fmt = 'o', ms=2.5 , color = self._filter_name.get(self.filter) )
        #     add_plot.plot(x_ , rise_typeII(x_, a_, t_exp_, n_), ls = '--' , label= n_string , color = _fitcolor )

        #     add_plot.plot(x_[0] , rise_typeII(x_, a_, t_exp_, n_)[0], color = 'white',label= chi_string )
        #     add_plot.plot(x_[1] , rise_typeII(x_, a_, t_exp_, n_)[1], color = 'white',label= chi_fit_string )


        #     add_plot.axvline(t_exp_, color = _fitcolor ,ls=':',alpha = 0.5, label = tx_string )
        #     # plt.axvspan(t_exp_ - dt_exp_, t_exp_ + dt_exp_, color = 'lightblue', alpha = 0.5)


        #     plt.xlabel('Time from ASFD [days]', size = 15)
        #     plt.ylabel('Flux [Jy]', size = 15)

        #     plt.legend(fontsize = 10)

        
        