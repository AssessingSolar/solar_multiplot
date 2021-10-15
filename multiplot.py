import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pvlib
import matplotlib.dates as mdates

comparison_flags = ['flag3lowSZA', 'flag3highSZA']
limit_flags = ['flagERLGHI','flagERLDIF','flagERLDNI','flagPPLGHI','flagPPLDIF','flagPPLDNI']
automatic_flags = ['flagKnKt','flagKn','flagKt','flagKlowSZA','flagKhighSZA','flagKKt','flag3lowSZA','flag3highSZA',
         'flagERLGHI','flagERLDIF','flagERLDNI','flagPPLGHI','flagPPLDIF','flagPPLDNI','flagTracker']

# Function for retrieving SRTM horzion profile (writte by Yves-Marie Saint-Drenan)
def wps_Horizon_SRTM(location):
    # location: geopoint [lat lon elev]

    import uuid
    import os
    import pandas
    from urllib.request import urlopen

    if (len(location) == 2):
        location[3] = -999;

    #unique id
    uid = uuid.uuid4().hex

    fic_output_csv = 'horizon_srtm_output_{}.csv'.format(uid)

    str_wps = 'http://toolbox.webservice-energy.org/service/wps?service=WPS&request=Execute&identifier=compute_horizon_srtm&version=1.0.0&DataInputs=';
    datainputs_wps = 'latitude={:.6f};longitude={:.6f};altitude={:.1f}'\
	.format(location[0], location[1], location[2]);
    
    #print(datainputs_wps)
    response = urlopen('{}{}'.format(str_wps,datainputs_wps))
    HZ = pandas.read_csv(response,delimiter=';',comment='#',header=None,skiprows=16,nrows=360,names=['AZIMUT', 'ELEVATION'])
     
    return HZ

def multi_plot(df, station_name, location, local_tz=None, kind='raw', save=False):
    df = df.copy()
    
    if local_tz != None:
        df.index = df.index.tz_convert(local_tz)
        
    if kind == 'final': # select only non-flagged and daytime data for final data
        df.loc[(df['SZA']<90)&(df['flagTotal']==1)|(df['flagManual']==1), ['GHI','DIF','DNI']] = np.nan

    if df.loc[df['ELV']>0, ['GHI','DNI','DIF']].dropna().empty:
        return

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    plt.rcParams['figure.constrained_layout.use'] = True

    layout = [['DNI_diff', 'GHI_hr_day', 'GHI_hr_day'],
              ['DNI_diff', 'DNI_hr_day', 'DNI_hr_day'],
              ['GHI_ratio', 'clearsky', 'clearsky'],
              ['GHI_scatter', 'K_Kt', 'GHI_elv_azi'],
              ['GHI_am_pm', 'Kn_Kt', 'DNI_elv_azi']]

    grid_ratios = {'height_ratios':[1,1,1,1.8,1.8], 'width_ratios':[1.3,1,2.2]}


    ################## HORIZON LINE ###########################################
    horizon_line = wps_Horizon_SRTM([location.latitude, location.longitude, location.altitude])


    ############### CORECTIONS FOR SOUTHERN HEMISPHERE #######################
    df['AZI_corr'] = df['AZI']
    if location.latitude<0:
        df['AZI_corr'] = (df['AZI'].add(180).mod(360))
        horizon_line['AZIMUT'] = horizon_line['AZIMUT'].add(180).mod(360)
        horizon_line = horizon_line.sort_values(by='AZIMUT', axis='rows')



    fig = plt.figure(constrained_layout=True, figsize=(20*1.3,15*1.3))
    axes = fig.subplot_mosaic(layout, gridspec_kw=grid_ratios)
    
    ########################### DNI DIFFERENCE PLOT ################################
    dni_diff_min = -70 # minimum ylimit
    dni_diff_max = 70 
    df_2d = df.groupby([df['AZI_corr'].round(0), df['DNI_diff'].round(0)])['DNI'].count().unstack('AZI_corr')
    df_2d = df_2d.reindex(np.arange(dni_diff_min, dni_diff_max), fill_value=np.nan) # Enforces limits and adds missing values
    df_2d = df_2d.divide(df_2d.sum()/100, axis='columns')

    extent = [df_2d.columns.min(), df_2d.columns.max(), df_2d.index.min(), df_2d.index.max()]
    im = axes['DNI_diff'].imshow(df_2d, aspect='auto', origin='lower',
                                 cmap='jet', extent=extent, vmin=0, vmax=df_2d.quantile(0.99).median())
    axes['DNI_diff'].set_ylabel('DNI$_{meas}$-DNI$_{calc}$ [W/m$^2$]')
    axes['DNI_diff'].set_xlabel('Azimuth [°N]')

    axes['DNI_diff'].grid(alpha=0.5)
    axes['DNI_diff'].set_ylim(dni_diff_min, dni_diff_max)
    if location.latitude < 0:
        axes['DNI_diff'].set_xticklabels(np.mod(axes['DNI_diff'].get_xticks()-180,360).astype(int))
    #cbar = plt.colorbar(im, ax=axes['DNI_diff'], orientation='vertical', pad=0.01,
    #                    label='Frequency within x bin [%]', ticks=np.arange(0,df_2d.max().max(),2))


    ################# SUNRISE SUNSET TIMES ##################
    days = pd.date_range(df.index[0], df.index[-1], tz=local_tz) # List of days for which to calculate sunrise/sunset

    sunrise_sunset = location.get_sun_rise_set_transit(days)

    # Convert sunrise/sunset from Datetime to hours (decimal)
    sunrise_sunset['sunrise'] = sunrise_sunset['sunrise'].dt.hour + sunrise_sunset['sunrise'].dt.minute/60
    sunrise_sunset['sunset'] = sunrise_sunset['sunset'].dt.hour + sunrise_sunset['sunset'].dt.minute/60


    #################### 2D plot of GHI & DNI where y-axis is hours and x-axis is day ####################
    #df_2d = df[['GHI','DNI']].set_index([df.index.date, df.index.hour+df.index.minute/60]).unstack(level=0)
    df_2d = df[['GHI','DNI']].set_index([df.index.date, df.index.hour*60+df.index.minute]).unstack(level=0)
    df_2d = df_2d.reindex(np.arange(1440)) # ensure all minutes are present
    df_2d.index = df_2d.index/60 # necessary to first convert afterwards, due to rounding issues.
    
    # Calculate the extents of the 2D plot, in the format [x_start, x_end, y_start, y_end]
    xlims = mdates.date2num([df.index[0].date(), df.index[-1].date()])
    extent = [xlims[0], xlims[1], 0, 24]
    
    for c in ['GHI_hr_day', 'DNI_hr_day']: # names of plot
        component = c.split('_')[0]
        im = axes[c].imshow(df_2d[component], 
                            aspect='auto', origin='lower', cmap='jet',
                            extent=extent, vmin=0, vmax=1100)#df[component].quantile(0.999))

        axes[c].set_xlim(xlims)
        axes[c].set_yticks(np.arange(0,25,6))
        axes[c].set_ylabel('Time of day [h]')
        axes[c].set_facecolor('grey')
        axes[c].xaxis_date()
        axes[c].xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))
        axes[c].plot(mdates.date2num(sunrise_sunset.index), sunrise_sunset[['sunrise', 'sunset']], 'm--')
        #axes[c].set_xlabel('Date')
        divider = make_axes_locatable(axes[c])
        cax1 = divider.append_axes('right', size='1%', pad=0.1)
        cbar1 = fig.colorbar(im, cax=cax1, orientation='vertical', pad=0.01, label='{} [W/m$^2$]'.format(component))


    ########################################## Clear sky plot  ###########################################
    df.loc[df['GHI']>50, 'cams_clearsky'].plot(ax=axes['clearsky'], style='bo', markersize=0.05, alpha=0.5, rot=0)
    axes['clearsky'].set_xlim([df.index[0], df.index[-1]]), axes['clearsky'].set_ylim([0,2])
    axes['clearsky'].xaxis.set_major_formatter(mdates.DateFormatter('%Y %b'))
    axes['clearsky'].set_ylabel('Clear sky index [-]')
    #axes['clearsky'].axhline(1, c='k', lw=0.5)
    divider2 = make_axes_locatable(axes['clearsky'])
    cax2 = divider2.append_axes("right", size="1%", pad=0.1)
    cax2.remove()

    ################################# GHI_meas / GHI_calc vs. time  ################################
    df['GHI'].divide(df['GHI_calc'])[(df['GHI']>50)&(df['SZA']<90)].plot(ax=axes['GHI_ratio'], style='o', c='grey',
                                        markersize=1, label='data GHI>50 W/m$^2$', rot=0)
    if kind=='automatic_flags':
        df['GHI'].divide(df['GHI_calc'])[(df['flagAutomatic']==1)&(df['SZA']<90)].plot(ax=axes['GHI_ratio'],
                                        style='o', c='red', markersize=1, label='flagged', rot=0)
    axes['GHI_ratio'].plot([0,1400], [0,1400], c='k', lw=0.5)
    axes['GHI_ratio'].set_ylim([0,4]), axes['GHI_ratio'].set_ylabel('GHI/GHI$_{calc}$')
    axes['GHI_ratio'].set_xlabel(''), axes['GHI_ratio'].set_xlim(df.index[0], df.index[-1])
    axes['GHI_ratio'].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    


    ################################ GHI_calc vs. GHI_meas ################################
    df[df['GHI']>50].plot.scatter(ax=axes['GHI_scatter'], x='GHI', y='GHI_calc', c='grey', s=1, label='data GHI>50 W/m$^2$')
    if kind=='automatic_flags':
        df[df[comparison_flags].max(axis='columns')==1].plot.scatter(ax=axes['GHI_scatter'], x='GHI', y='GHI_calc', c='red', s=1, label='flagged')
    axes['GHI_scatter'].set_xlim(0,1400), axes['GHI_scatter'].set_ylim(0, 1400)
    axes['GHI_scatter'].set_xlabel('GHI$_{meas}$'), axes['GHI_scatter'].set_ylabel('GHI$_{calc}$')
    axes['GHI_scatter'].get_legend().remove()


    ################################ Scatter plots for K vs Kt plot ################################
    a = 1
    if kind=='final':
        df[(df['flagTotal']!=1)&(df['SZA']<90)].plot.scatter(ax=axes['K_Kt'], x='Kt', y='K', c='grey', s=0.5, alpha=a, label='data')
        df[(df['flagTotal']!=1)&(df['SZA']<90)&(df['GHI']<=50)].plot.scatter(ax=axes['K_Kt'], x='Kt', y='K', c='black', s=0.5, alpha=a, label='data ghi below 50 W/m$^2$')
    elif kind=='automatic_flags':
        df[(df['GHI']>50)&(df['SZA']<90)].plot.scatter(ax=axes['K_Kt'], x='Kt', y='K', c='grey', s=0.5, alpha=a, label='data')
        df[(df['GHI']<=50)&(df['SZA']<90)].plot.scatter(ax=axes['K_Kt'], x='Kt', y='K', c='black', s=0.5, alpha=a, label='data ghi below 50 W/m$^2$')
        df[df['k_flags']==1].plot.scatter(ax=axes['K_Kt'], x='Kt', y='K', c='red', s=1, alpha=a, label='flagged')
    else:
        df[(df['SZA']<90)].plot.scatter(ax=axes['K_Kt'], x='Kt', y='K', c='grey', s=0.5, alpha=a, label='data ghi below 50 W/m$^2$')
    
    axes['K_Kt'].plot([0,0.6,0.6,1.35,1.35], [1.1,1.1,0.96,0.96,0], c='blue', label='limits')
    axes['K_Kt'].plot([0,0.6], [1.05,1.05], c='blue', linestyle='dashed', label='limits')
    axes['K_Kt'].set_xlim(0, 1.4), axes['K_Kt'].set_ylim(0,1.2)
    axes['K_Kt'].set_xlabel('Kt = GHI / ETH [-]'), axes['K_Kt'].set_ylabel('K = DHI / GHI [-]')
    axes['K_Kt'].get_legend().remove()


    ################################# Scatter plots for Kn vs Kt plot ################################
    if kind=='final': # plot only unflagged data
        df[(df['SZA']<90)&(df['flagTotal']!=1)].plot.scatter(ax=axes['Kn_Kt'], x='Kt', y='Kn', c='grey', s=0.5, alpha=a, label='data')
    elif kind=='automatic_flags': # plot all data in grey and flagged in red
        df.plot.scatter(ax=axes['Kn_Kt'], x='Kt', y='Kn', c='grey', s=0.5, alpha=a, label='data')
        df[df['k_flags']==1].plot.scatter(ax=axes['Kn_Kt'], x='Kt', y='Kn', c='red', s=1, alpha=a, label='flagged')
    else: # plot all data
        df[df['SZA']<90].plot.scatter(ax=axes['Kn_Kt'], x='Kt', y='Kn', c='grey', s=0.5, alpha=a, label='data')

    axes['Kn_Kt'].plot([0,0.96,1.35,1.35], [0,0.95,0.95,0], c='blue', label='limits')
    axes['Kn_Kt'].set_xlim(0, 1.4), axes['Kn_Kt'].set_ylim(0,1)
    axes['Kn_Kt'].set_xlabel('Kt = GHI / ETH [-]'), axes['Kn_Kt'].set_ylabel('Kn = DNI / ETN [-]')
    axes['Kn_Kt'].get_legend().remove()
    

    ####### Maximum GHI/DNI values (solar elevation vs. azimuth) ######### very slow

    df_2d = df.loc[df['ELV']>=0, ['GHI','DNI']].groupby([df['AZI_corr'].round(0), df['ELV'].multiply(2).round(0).divide(2)]).max().unstack(level=0).sort_index()
    #extent = [0, 360, 0, df_2d.index.max()] # This gave incorrect x-axis as not all azimuth values are present
    for c in ['GHI_elv_azi','DNI_elv_azi']:
        component = c.split('_')[0]
        extent = [df_2d[component].columns.min(), df_2d[component].columns.max(), df_2d[component].index.min(), df_2d[component].index.max()]
        im = axes[c].imshow(df_2d.loc[0:, component], aspect='auto', origin='lower', cmap='jet',
                            extent=extent, vmin=0, vmax=df[component].quantile(0.99))

        
        axes[c].grid(alpha=0.5)
        axes[c].set_xlim(0,360)
        axes[c].set_xticks(np.arange(0,360+60,60))
        axes[c].plot(horizon_line['AZIMUT'], horizon_line['ELEVATION'], c='k', lw=0.5)
        if location.latitude<0:
            axes[c].set_xticklabels(np.mod(axes[c].get_xticks()+180,360))

        axes[c].set_yticks(np.arange(0,round(df_2d.index.max()/10)*10+10+10,10))
        axes[c].set_ylim(0,None)
        axes[c].set_ylabel('Solar elevation [°]')
        axes[c].set_xlabel('Solar azimuth [°N]')
        divider3 = make_axes_locatable(axes[c])
        cax3 = divider3.append_axes("right", size="2%", pad=0.1)
        cbar3 = fig.colorbar(im, cax=cax3, pad=0.01, orientation='vertical', label='Max. {} [W/m$^2$]'.format(component))


    ####################### GHI AM/PM ratio plot #######################
    df['AZI_SOUTH'] = (df['AZI']-180).abs().divide(2).round(0).multiply(2)
    df['GHI_am'] = df.loc[df.index.hour>=12, 'GHI']
    df['GHI_pm'] = df.loc[df.index.hour<=11, 'GHI']

    # Calculate ratio of GHI before and outer AZI=SOUTH for every day
    df_am_pm = df.loc[(df['ELV']>1) & (df['GHI']>50)].groupby(['date','AZI_SOUTH']).mean()[['GHI_am','GHI_pm']].reset_index()
    df_am_pm['GHI_am_pm'] = df_am_pm['GHI_am'].divide(df_am_pm['GHI_pm'], axis='rows').multiply(200).round(0).divide(200)
    df_am_pm.loc[df_am_pm['GHI_am_pm'].abs()-1>0.1, 'GHI_am_pm'] = np.nan
    df_am_pm = df_am_pm.reset_index().dropna().groupby(['AZI_SOUTH','GHI_am_pm']).count()
    df_am_pm = df_am_pm['date'].unstack('AZI_SOUTH')[0.95:1.05]
    df_am_pm = df_am_pm.divide(df_am_pm.sum().sum())*100

    extent = [df_am_pm.columns.min(), df_am_pm.columns.max()+1, df_am_pm.index.min(), df_am_pm.index.max()]
    im = axes['GHI_am_pm'].imshow(df_am_pm, aspect='auto', extent=extent, vmax=None, cmap='plasma')
    axes['GHI_am_pm'].set_xlim(0, df_am_pm.columns.max()+1)
    axes['GHI_am_pm'].set_xlabel('Abs. azimuth [°S]'), axes['GHI_am_pm'].set_ylabel('GHI$_{AM}$ / GHI$_{PM}$')

    # Set title and meta-data
    ghi_kwh = df['GHI'].resample('1h').mean().clip(lower=0).sum()/1000 # GHI sum kwh/m2
    dni_kwh = df['DNI'].resample('1h').mean().clip(lower=0).sum()/1000 # DNI sum kwh/m2
    fig.suptitle('Station name: {:12}  Latitude: {:06.3f}°N  Longitude: {:>4.3f}°E  First date: {}  Last date: {}  GHI sum={:>4.0f} kWh/m$^2$  DNI sum={:04.0f} kWh/m$^2$'.format(
                 station_name, location.latitude, location.longitude, df.index[0].date(), df.index[-1].date(), ghi_kwh, dni_kwh), fontsize=16, y=0.91)
    
    if save==True:
        fig.savefig('Output_plots/{}_{}.png'.format(station_name, kind), dpi=250, bbox_inches='tight')
        plt.close()




