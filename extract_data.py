
Hello there

Hi German. From Rob. 11:32 am

This is rob's second change. 11:42am'

"Hello everyone, from Dan, test add"

def useful_PNiMo(my_exp, pNiMos, all_flux):
    useful_pNiMo = 0
    for pNiMo, flux in zip(pNiMos, all_flux):
        if pNiMo > 0.75 and pNiMo < 2 and flux > 1e8:
            useful_pNiMo += 1
    useful_pNiMo = 100 * (useful_pNiMo / float(my_exp.motor_driver.number_of_bursts))
    return useful_pNiMo


def slope_of_PNiMos(pNiMos, all_flux):
    all_values = []
    x = []
    i = 0
    for pNiMo, flux in zip(pNiMos, all_flux):
        if pNiMo > 0.75 and pNiMo < 2 and flux > 1e8:
            all_values.append(pNiMo)
            x.append(i)
        i += 1
    if len(x) > 0:
        z = np.polyfit(x, all_values, 1)
        slope = z[0]
    else:
        slope = float('nan')
    return slope
    
def std_pNiMo(pNiMos, all_flux):
    all_values = []
    x = []
    i = 0
    for pNiMo, flux in zip(pNiMos, all_flux):
        if pNiMo > 0.75 and pNiMo < 2 and flux > 1e8:
            all_values.append(pNiMo)
            x.append(i)
        i += 1
    if len(x) > 0:
        z = np.polyfit(x, all_values, 1)
        p = np.poly1d(z)
        residuals = abs(all_values - p(x))
        std = np.std(residuals)
    else:
        std = float('nan')
    return std
    
    
Tk().withdraw()
path = askopenfilename(title='Pick Analysis File or Zipfile of the Experiment')
print path
save_path = askopenfilename(title='Pick Save Location')
my_exp = Experiment(path, pseudo_NiMo_step_size=1)
pNiMos = my_exp.pseudo_NiMo_burst_evolution
all_flux = my_exp.average_flux_within_bursts


avg_flux = my_exp.average_flux_whole_run
avg_flux_above_threshold = my_exp.average_flux_above_flux_threshold
useful_pnimo = useful_PNiMo(my_exp, pNiMos, all_flux)
slope_pnimo = slope_of_PNiMos(pNiMos, all_flux)
std_pnimo = std_pNiMo(pNiMos, all_flux)
nbr_of_bursts = my_exp.motor_driver.number_of_bursts
nbr_of_bursts_above_1e8 = my_exp.bursts_over_flux_threshold
#counts, w_intensity, R_2 = sct.calculate_W_intensity(my_exp.x_ray_detector.spectrum, my_exp.x_ray_detector.spectrum_axis)
exp_name = my_exp.test_parameters.test_number

if "Cu_filter" in my_exp.test_parameters.test_number:
    peak_energy = my_exp.x_ray_detector.peak_energy_cu_filter
else:
    peak_energy = float('nan')

data = {' Experiment': exp_name,
           'Avg Flux whole run': [avg_flux],
           'Avg Flux Above Threshold': [avg_flux_above_threshold],
           'Useful PNiMo': [useful_pnimo],
           'Slope PNiMo': [slope_pnimo],
           'Std PNiMo': [std_pnimo],
           'Nbr of Bursts': [nbr_of_bursts],
           'Peak Energy': [peak_energy], 
           'Flux above 1e8': [nbr_of_bursts_above_1e8],
           }#'W-intensity': [w_intensity] }



if not save_path == '':
    df = pd.read_csv(save_path, sep='\t')
    dic = df.to_dict(orient='list')
    for key in dic.keys():
        dic[key].append(data[key][0])
else:
    dic = data
    save_path = 'results' + '.csv'
    
df = pd.DataFrame.from_dict(dic)
df.to_csv(save_path, sep= '\t', index=False)






