#!/usr/bin/env python
# coding: utf-8

# In[1]:


import proposal as pp
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from pandas.core.common import flatten
import vg


# In[2]:


pp.InterpolationDef.path_to_tables = "$HOME/.cache/PROPOSAL"


# # Define functions

# In[3]:


def propagate_deflected_muons(initial_energies, minimum_energies, inter_type, deflection_type='continuous+stochastic', e_cut=500,                              v_cut=0.05, cont_rand=False, scattering_method="highlandintegral", beta=1.0, rnd_seed=1337, initial_direction=[0, 0, 1]):

    '''Propagate muon tracks with deflection. Scaling of Bremsstrahlung opening angle can be done by beta.
    
    Parameters
    ----------
    initial_energies: list of energies
    minimum_energs: list of energies, lower propagation limit
    inter_type: list of interaction types for propagation/deflection
    deflection_type: string, choose one:
            1. 'contiuous+stochastic'
            2. 'continuous'
            3. 'stochastic'
    beta: scaling factor for Bremsstrahlung
    e_cut, v_cut, cont_rand: usual PROPOSAL energy cut settings
    initial_direction: list of initial direction (cartesian coordinates)
    '''
    pp.RandomGenerator.get().set_seed(rnd_seed)
    args = {
            "particle_def": pp.particle.MuMinusDef(),
            "target": pp.medium.Ice(),
            "interpolate": True,
            "cuts": pp.EnergyCutSettings(e_cut, v_cut, cont_rand)
            }

    cross = pp.crosssection.make_std_crosssection(**args)
    multiple_scatter = pp.make_multiple_scattering(scattering_method, args["particle_def"], args["target"], cross, True)
    stochastic_deflect = pp.make_default_stochastic_deflection(inter_type,
            args["particle_def"], args["target"])


    collection = pp.PropagationUtilityCollection()
    collection.displacement = pp.make_displacement(cross, True)
    collection.interaction = pp.make_interaction(cross, True)
    collection.time = pp.make_time(cross, args["particle_def"], True)
    collection.decay = pp.make_decay(cross, args["particle_def"], True)

    if deflection_type == 'stochastic':
        print('stochastic deflection')
        if pp.particle.Interaction_Type.brems in inter_type:
            collection.scattering = pp.scattering.ScatteringMultiplier(stochastic_deflect, [(pp.particle.Interaction_Type.brems, beta)])
        else:
            collection.scattering = pp.scattering.ScatteringMultiplier(stochastic_deflect, [(inter_type[0], 1.0)])
    elif deflection_type == 'continuous':
        print('continuous deflection')
        collection.scattering = pp.scattering.ScatteringMultiplier(multiple_scatter, 1.0)
    elif deflection_type == 'continuous+stochastic':
        print('continuous and stochastic deflection')
        if pp.particle.Interaction_Type.brems in inter_type:
            collection.scattering = pp.scattering.ScatteringMultiplier(multiple_scatter, stochastic_deflect, 1.0, [(pp.particle.Interaction_Type.brems, beta)])
        else:
            collection.scattering = pp.scattering.ScatteringMultiplier(multiple_scatter, stochastic_deflect, 1.0, [(inter_type[0], 1.0)])

    utility = pp.PropagationUtility(collection = collection)
    detector = pp.geometry.Sphere(pp.Vector3D(0,0,0), 1e20)
    density_distr = pp.density_distribution.density_homogeneous(args["target"].mass_density)

    
    prop = pp.Propagator(args["particle_def"], [(detector, utility, density_distr)])

    init_state = pp.particle.ParticleState()
    init_state.position = pp.Vector3D(0, 0, 0)
    init_state.direction = pp.Vector3D(initial_direction[0], initial_direction[1], initial_direction[2])

    tracks = []
    for E_i, E_min in zip(tqdm(initial_energies), minimum_energies):
        init_state.energy = E_i # initial energy in MeV
        track = prop.propagate(init_state, max_distance = 1e9, min_energy = E_min)
        tracks.append(track)
        
    return tracks
                

def get_opening_angles_and_types(tracks):
    '''Calculate opening angles of each interaction, maximum energy losses and its types.
    
    Parameters
    ----------
    tracks: list of propagated muon tracks
        
    Returns
    ------
    opening_type_dict: dictionary, see documentation in return below
    
    '''
    types = {
            'Interaction_Type.epair': 0,
            'Interaction_Type.brems': 0,
            'Interaction_Type.photonuclear': 0,
            'Interaction_Type.ioniz': 0,
    }
    stochastic_opening = []
    continuous_opening = []
    all_types = []
    types_no_cont = []
    outcoming_angle = []
    st_relative_loss = []
    st_loss = []
    co_loss = []
    all_relative_losses = []
    initial_energies = []
    for track in tqdm(tracks):  
        st_op_temp = []
        co_op_temp = []
        all_types_temp = []
        types_no_cont_temp = []
        st_relative_loss_temp = []
        st_loss_temp = []
        co_loss_temp = []
        all_rel_loss_temp = []
        theta_last = track.track_directions()[1].theta
        phi_last = track.track_directions()[1].phi
        E_i = track.track_energies()[0]
        E_last = track.track_energies()[1]
        for typ,angle,E in zip(track.track_types()[2:], track.track_directions()[2:],                                track.track_energies()[2:]):
            angle_op = np.rad2deg(get_angle_deviation(phi_last, theta_last, angle.phi,                                                 angle.theta))
            if str(typ) in types:
                types[str(typ)] += 1
                st_op_temp.append(angle_op)
                all_types_temp.append(str(typ))
                types_no_cont_temp.append(str(typ))
                st_relative_loss_temp.append((E_last - E) / E_last)
                st_loss_temp.append(E_last - E)
            else:
                co_op_temp.append(angle_op)
                all_types_temp.append(str(typ))
                co_loss_temp.append(E_last - E)
            theta_last = angle.theta
            phi_last = angle.phi
            all_rel_loss_temp.append((E_last - E)/E_i)
            E_last = E
            
        stochastic_opening.append(st_op_temp)
        continuous_opening.append(co_op_temp)
        all_types.append(all_types_temp)
        types_no_cont.append(types_no_cont_temp)
        angle_out = np.rad2deg(get_angle_deviation(phi_last, theta_last,                             track.track_directions()[1].phi, track.track_directions()[1].theta))
        outcoming_angle.append(angle_out)
        st_relative_loss.append(st_relative_loss_temp)
        st_loss.append(st_loss_temp)
        co_loss.append(co_loss_temp)
        all_relative_losses.append(all_rel_loss_temp)
        initial_energies.append(E_i)
        
    
    opening_typ_dict = {
        'stoch_opening': stochastic_opening, # in degree
        'cont_opening': continuous_opening, # in degree
        'types_counter': types, # counter of each interaction
        'all_types': all_types, # list of different interaction types
        'types_no_cont': types_no_cont, # list of different interactions types without cont losses
        'outcoming_angle': outcoming_angle, # opening angle between intial direction and end direction 
        'st_relative_loss': st_relative_loss, # relative energy losses of stochastic interactions (normed on current energy)
        'st_loss': st_loss, # stochastic energy losses
        'co_loss': co_loss, # continuous energy losses
        'all_relative_losses_in': all_relative_losses, # all relative energy losses (normed on initial energy)
        'initial_energies': initial_energies, # initial energies in MeV
    }
    return opening_typ_dict
    
       
def get_angle_deviation(azimuth1, zenith1, azimuth2, zenith2):
    """Get opening angle of two vectors defined by (azimuth, zenith)
    Parameters
    ----------
    azimuth1 : 
        Azimuth of vector 1 in rad.
    zenith1 : 
        Zenith of vector 1 in rad.
    azimuth2 : 
        Azimuth of vector 2 in rad.
    zenith2 : 
        Zenith of vector 2 in rad.
    Returns
    -------
        The opening angle in rad between the vector 1 and 2.
        Same shape as input vectors.
    """
    cos_dist = (np.cos(azimuth1 - azimuth2) *
                np.sin(zenith1) * np.sin(zenith2) +
                np.cos(zenith1) * np.cos(zenith2))
    cos_dist = np.clip(cos_dist, -1., 1.)
    return np.arccos(cos_dist)  

def get_maximum_loss_types_and_deflections(st_loss, types_no_cont, stoch_opening, E_i, co_loss, cont_opening):
    '''Get maximum interaction types and its corresponding deflection angles.
    
    Parameters
    ----------
    st_loss: array_like, stochastic losses
    types_no_cont: array_like, stochastic interaction types
    stoch_opening: array_like, stochastic deflections angles
    E_i: array_like, initial muon energies
    co_loss: array_like, continuous energy losses
    cont_opening: array_like, continuous deflection angles
    
    Returns
    -------
    max_dict: dictionary, see below
    '''
    max_loss_types = []
    max_defl = []
    max_rel_loss_in = []
    max_loss = []
    co_max_defl = []
    co_max_rel_loss_in = []
    co_max_loss = []
    
    for i in range(len(st_loss)):
        if len(st_loss[i]) > 0:
            arg_max = np.argmax(st_loss[i])
            max_loss_types.append(types_no_cont[i][arg_max])
            max_defl.append(stoch_opening[i][arg_max])
            max_rel_loss_in.append(np.max(st_loss[i])/E_i[i])
            max_loss.append(np.max(st_loss[i]))
        else:
            max_defl.append(0)
            max_rel_loss_in.append(0)
            max_loss.append(0)
            
    
    for i in range(len(co_loss)):
        if len(co_loss[i]) > 0:
            arg_max = np.argmax(co_loss[i])
            co_max_defl.append(cont_opening[i][arg_max])
            co_max_rel_loss_in.append(np.max(co_loss[i])/E_i[i])
            co_max_loss.append(np.max(co_loss[i]))
        else:
            co_max_defl.append(0)
            co_max_rel_loss_in.append(0)
            co_max_loss.append(0)
        
    max_dict = {
        'max_loss_types': max_loss_types, # type of maximum energy loss
        'max_defl': max_defl, # deflection of maximum stochastic loss
        'max_rel_loss_in': max_rel_loss_in, # relative maximum energy loss normed on initial energy
        'max_loss': max_loss, # absolute maximum stochastic energy loss
        'co_max_defl': co_max_defl, # deflection of maximum continuous energy loss
        'co_max_rel_loss_in': co_max_rel_loss_in, # maximum continuous relative loss normed on initial energy
        'co_max_loss': co_max_loss, # absolute maximum continuous energy loss
    }
    return max_dict


# In[4]:


###### settings #####
inter_type =[pp.particle.Interaction_Type.ioniz, pp.particle.Interaction_Type.brems,             pp.particle.Interaction_Type.photonuclear, pp.particle.Interaction_Type.epair]
v_cut = 0.05
scattering_method = "moliere"
beta = 1.0
deflection_type = 'continuous+stochastic'
initial_direction = [0, 0, 1]

# set energies
number_tracks = 100
rnd = np.random.RandomState(42)
E = 1e8
E_i = rnd.exponential(E, size=number_tracks) 
E_min = np.ones(number_tracks) * 5e5 # set lowest energy to 500 GeV
#####################

tracks = propagate_deflected_muons(initial_energies=E_i, minimum_energies=E_min, inter_type=inter_type,             deflection_type=deflection_type, v_cut=v_cut, scattering_method=scattering_method,             beta=beta, initial_direction=initial_direction)


# In[5]:


bins = np.logspace(-3, 8, 30)
plt.hist(np.array(E_i)/1e3, bins=bins, histtype='step', label='initial energy')
plt.hist(E_min/1e3, bins=bins, histtype='step', label='outcoming energy')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$E\,/\,$GeV')
plt.ylabel('counts')
plt.legend()


# # Work with opening angles

# In[6]:


opening_typ_dict = get_opening_angles_and_types(tracks)


# In[7]:


normed_types_of_all_interactions = {k: np.round(opening_typ_dict['types_counter'][k]/sum(opening_typ_dict['types_counter'].values()), 2) for k in opening_typ_dict['types_counter']}
normed_types_of_all_interactions


# In[8]:


max_types = {
            'Interaction_Type.epair': 0,
            'Interaction_Type.brems': 0,
            'Interaction_Type.photonuclear': 0,
            'Interaction_Type.ioniz': 0,
            'Interaction_Type.continuousenergyloss': 0,
    }
max_types_no_cont = {
            'Interaction_Type.epair': 0,
            'Interaction_Type.brems': 0,
            'Interaction_Type.photonuclear': 0,
            'Interaction_Type.ioniz': 0,
    }

for loss,typ in zip(opening_typ_dict['all_relative_losses_in'], opening_typ_dict['all_types']):
    if len(loss) < 1:
        continue
    max_types[str(typ[np.argmax(loss)])] += 1
    no_cont_typ = [x for x in typ if x != 'Interaction_Type.continuousenergyloss']
    no_cont_loss = [y for x,y in zip(typ,loss) if x != 'Interaction_Type.continuousenergyloss']
    max_types_no_cont[str(no_cont_typ[np.argmax(no_cont_loss)])] += 1


print('all maximum loss types (stoch+cont): \n', max_types, '\n')
normed_max_types = {k: np.round(max_types[k]/sum(max_types.values()), 2) for k in max_types}
print('normed all maximum loss types (stoch+cont): \n', normed_max_types, '\n')

print('all maximum loss types (stoch): \n', max_types_no_cont, '\n')
normed_max_types_no_cont = {k: np.round(max_types_no_cont[k]/sum(max_types_no_cont.values()), 2) for k in max_types_no_cont}
print('normed all maximum loss types (stoch): \n', normed_max_types_no_cont, '\n')


# In[9]:


bins = np.logspace(-9, 1, 30)
plt.hist(list(flatten(opening_typ_dict['stoch_opening'])), bins=bins, histtype='step', label='stoch op')
plt.hist(list(flatten(opening_typ_dict['cont_opening'])), bins=bins, histtype='step', label='cont op')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('opening angle of each interaction in degree')
plt.ylabel('counts')
plt.legend()
events_s = len([i for i in opening_typ_dict['stoch_opening'] if len(i)>0])
events_c = len([i for i in opening_typ_dict['cont_opening'] if len(i)>0])
plt.title('stoch events: {}, cont events: {}'.format(events_s, events_c))



print('minimum values for opening calculation:')
nu = [0.999999, 0.999999, 0.9999999, 0.99999999]
test_dev = [get_angle_deviation(0.8*n, 0.7, 0.8, 0.7) for n in nu]
print(np.min(np.rad2deg(test_dev)))
print(np.rad2deg(test_dev))
print('---> numeric limit at ~1e-6')


# In[10]:


bins = np.logspace(-4, 1, 12)
n, x, _ = plt.hist(list(flatten(opening_typ_dict['outcoming_angle'])), bins=bins, histtype='step')
plt.xscale('log')
plt.xlabel('outcoming angle in degree')
plt.ylabel('counts')
plt.title('number events: {}'.format(np.sum(n)))


# In[11]:


bins = np.logspace(-9, 1, 30)
for typ in tqdm(['Interaction_Type.epair','Interaction_Type.ioniz',            'Interaction_Type.photonuclear','Interaction_Type.brems']):
    stoch_defl = []
    for j in range(len(opening_typ_dict['stoch_opening'])):
        helper = opening_typ_dict['stoch_opening'][j]
        st_defl = [helper[i] for i,t in enumerate(opening_typ_dict['types_no_cont'][j]) if t==typ]
        stoch_defl.append(st_defl)
    plt.hist(list(flatten(stoch_defl)), bins=bins, histtype='step', label=typ)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('opening angle of each interaction in degree')
plt.ylabel('counts')
plt.legend(bbox_to_anchor=(1,1))


# In[12]:


bins = (np.logspace(-7.5, 1, 31), np.logspace(-7, 0, 30))

x = list(flatten(opening_typ_dict['stoch_opening']))
y = list(flatten(opening_typ_dict['st_relative_loss']))

plt.hist2d(x, y, bins=bins, cmin=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta\,/\,°$')
plt.ylabel('relative energy loss') # normed on incoming energy
plt.title('all stochastic deflections')
plt.colorbar()


# In[13]:


max_dict = get_maximum_loss_types_and_deflections(opening_typ_dict['st_loss'],                                     opening_typ_dict['types_no_cont'], opening_typ_dict['stoch_opening'],                                     opening_typ_dict['initial_energies'], opening_typ_dict['co_loss'],                                     opening_typ_dict['cont_opening'])

for t in np.unique(max_dict['max_loss_types']):
    print(t,': ', max_dict['max_loss_types'].count(t), '({})'          .format(np.round(max_dict['max_loss_types'].count(t)/len(max_dict['max_loss_types']), 2)))


# In[14]:


bins = (np.logspace(-7, 2, 30), np.logspace(-5, 0, 30))
plt.hist2d(max_dict['max_defl'], max_dict['max_rel_loss_in'], bins=bins, cmin=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta\,/\,°$')
plt.ylabel('maximum relative energy loss') # normed on incoming energy
plt.title('stochastic losses')
plt.colorbar()


# In[15]:


bins = (np.logspace(-7, 2, 30), np.logspace(-5, 0, 30))
plt.hist2d(max_dict['co_max_defl'], max_dict['co_max_rel_loss_in'], bins=bins, cmin=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta\,/\,°$')
plt.ylabel('maximum relative energy loss') # normed on incoming energy
plt.title('continuous losses')
plt.colorbar()


# In[16]:


bins = (np.logspace(2, 10, 30), np.logspace(3, 5, 30))
plt.hist2d(max_dict['max_loss'], max_dict['co_max_loss'], bins=bins, cmin=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('energy of maximum loss (stochastic) in MeV')
plt.ylabel('energy of maximum loss (continuous) in MeV') # normed on incoming energy
plt.title('absolute energies')
plt.colorbar()


# In[17]:


phis = []
phi_diff = []
for track in tqdm(tracks):
    phi_last = track.track_directions()[1].phi
    phi_d = []
    for d in track.track_directions()[2:]:
        phi_d.append(d.phi - phi_last)
        phi_last = d.phi
    phi_diff.append(phi_d)
    phi = [d.phi for d in track.track_directions()]
    phis.append(phi)


# In[18]:


bins = np.linspace(-np.pi, np.pi, 21)
plt.hist(list(flatten(phis)), bins=bins, histtype='step', label='cont+stoch')
plt.xlabel('phi in rad')
plt.ylabel('counts')
plt.legend()


# In[19]:


bins = np.linspace(-2*np.pi, 2*np.pi, 21)
plt.hist(list(flatten(phi_diff)), bins=bins, histtype='step', label='cont+stoch')
plt.xlabel('phi_old - phi_new in rad')
plt.ylabel('counts')
plt.yscale('log')
plt.legend()


# In[ ]:




