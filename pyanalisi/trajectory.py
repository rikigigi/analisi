import numpy as np
import .pyanalisi as pa
import time

class Trajectory:
    KEYS = {'positions', 'cells'}
    KEYS_OPT = {'velocities', 'energy_constant_motion', 'ionic_temperature',
                'pressure', 'electronic_kinetic_energy', 'forces', 'stresses'}
    TIME_K = 'times'
    STEP_K = 'steps'

    DEFAULT_SAVE_REFERENCE = True

    @classmethod
    def from_lammps_binary(cls, fname, dt=None, traj_skip=1, symbols=None, wrap=False, nsteps=0, pk=None):
        '''
        Read a binary lammps file. traj_skip loads a step every traj_skip steps.
        Note that this is done at the numpy python level, and the trajectory object still has all steps inside it
        '''
        t = pa.Traj(fname)
        t.setWrapPbc(wrap)
        step_max=t.getNtimesteps()
        if nsteps > 0:
            if start+nstep > step_max:
                nsteps = step_max-start
        else:
            nsteps = step_max
        t.setAccessWindowSize(nsteps)
        t.setAccessStart(start)

        return cls(t,symbols_=symbols,dt=dt, pk=pk,traj_skip=traj_skip)

    @staticmethod
    def _copy(d, skip):
        return np.copy(d[::skip], order='C') if skip > 1 else d

    @staticmethod
    def analisi_cell2box(box_lammps):
        '''returns cell vectors arranged in ROWS.
           Analisi code always rotate the system to have a triangular cell matrix.'''
        cells = np.zeros((box_lammps.shape[0], 3, 3))
        for i in range(box_lammps.shape[0]):
            lc = box_lammps[i]
            if lc.shape[0] == 9:  # note that this matrix is transposed: the cell vectors are arranged in ROWS
                #                    a_x 0   0
                #                    b_x b_y 0
                #                    c_x c_y c_z
                cells[i] = np.array([[lc[3] * 2, 0, 0],
                                     [lc[6], lc[4] * 2, 0],
                                     [lc[7], lc[8], lc[5] * 2]
                                     ])
            elif lc.shape[
                0] == 6:  # for orthorombic cells there is no difference between row arranged and column arranged cell vectors
                cells[i] = np.array([[lc[3] * 2, 0, 0],
                                     [0, lc[4] * 2, 0],
                                     [0, 0, lc[5] * 2]])
            else:
                raise IndexError('wrong shape of cell array')
        return cells

    def _init_from_pa_traj(self, t, symbols_=None, traj_skip=1, save_reference=DEFAULT_SAVE_REFERENCE):
        if symbols_ is None:
            symbols_ = {}
        self.data['symbols'] = [symbols_[x] if x in symbols_ else str(x) for x in t.get_lammps_type().tolist()]
        self.symbols = self.data['symbols']
        self.data['positions'] = Trajectory._copy(t.get_positions_copy(), traj_skip)
        self.data['velocities'] = Trajectory._copy(t.get_velocities_copy(), traj_skip)
        box_lammps = Trajectory._copy(t.get_box_copy(), traj_skip)
        self.data['cells'] = Trajectory.analisi_cell2box(box_lammps)
        if save_reference:
            self.paTraj = t


    def __init__(self, t, symbols_=None, dt=None, pk=None, traj_skip=1, save_reference=DEFAULT_SAVE_REFERENCE):
        self.paTrajWrapped = None
        self.paTraj = None
        if symbols_ is None:
            symbols_ = {}

        if pk is None:
            pk = str(time.time_ns())

        self.data = {}
        if isinstance(t, list): ##list of whatever supported object
            nlist = [ Trajectory(it) for it in t ]
            self.symbols = nlist[0].symbols
            for key in self.KEYS:
                self.data[key] = Trajectory._copy(np.concatenate([it.get_array(key) for it in nlist], axis=0),
                                                  traj_skip)
            for key in self.KEYS_OPT:
                if key in t[0].get_arraynames():
                    self.data[key] = Trajectory._copy(np.concatenate([it.get_array(key) for it in nlist], axis=0),
                                                      traj_skip)
        elif isinstance(t, (pa.Traj, pa.Trajectory)):  ##pyanalisi trajectory
            self._init_from_pa_traj(t, symbols_=symbols_, traj_skip=traj_skip, save_reference=save_reference)
        elif isinstance(t, dict): ##simple dict with the data
            self.symbols = t['symbols']
            for key in self.KEYS:
                self.data[key] = Trajectory._copy(np.array(t[key]), traj_skip)
            for key in self.KEYS_OPT:
                if key in t:
                    self.data[key] = Trajectory._copy(np.array(t[key]), traj_skip)
        else: ##Trajectory-like interface, similar to aiida one
            try:
                self.symbols = t.symbols
                for key in self.KEYS:
                    self.data[key] = Trajectory._copy(t.get_array(key), traj_skip)
                for key in self.KEYS_OPT:
                    if key in t.get_arraynames():
                        self.data[key] = Trajectory._copy(t.get_array(key), traj_skip)
                if Trajectory.TIME_K in t.get_arraynames() and dt is None:
                    times=t.get_array(Trajectory.TIME_K)
                    dt=times[1]-times[0]
            except:
                raise RuntimeError(f'first argument cannot be {str(t)}')
        self.data[Trajectory,STEP_K] = np.arange(0, self.data['positions'].shape[0])
        if dt is None:
            dt=1.0
        self.data[Trajectory.TIME_K] = self.data[Trajectory.STEP_K] * dt * traj_skip
        self.dt = dt * traj_skip
        self.numsites = self.data['positions'].shape[1]
        self.pk = pk
        self.numsteps = self.data['positions'].shape[0]


    def get_array(self, name):
        if name in self.data:
            return self.data[name]
        else:
            raise KeyError(f'cannot find key {name}')

    def get_attribute(self, k):
        if k == 'symbols':
            return self.symbols
        else:
            raise KeyError(f'attribute {k} not present')

    def get_arraynames(self):
        return self.data.keys()

    @staticmethod
    def get_types_id_array(types_array):
        types = {}
        typeid = 0
        res = []
        for t in types_array:
            if not t in types:
                types[t] = typeid
                typeid = typeid + 1
            res.append(types[t])
        return np.array(res, dtype='int32')

    def get_analisi_traj(self, tskip=1, wrapped=False, save_reference=DEFAULT_SAVE_REFERENCE):
        if self.paTraj and tskip==1 and not wrapped:
            return self.paTraj
        elif self.paTrajWrapped and tskip==1 and  wrapped:
            return self.paTrajWrapped

        if tskip < 1:
            raise IndexError(f'cannot have a trajectory skip of {tskip}')
        pos = self.get_array('positions')[::tskip].copy(order='C')
        if 'velocities' in self.get_arraynames():
            vel = self.get_array('velocities')[::tskip].copy(order='C')
        else:
            vel = np.zeros(pos.shape)
        cel = self.get_array('cells')[::tskip].transpose((0, 2, 1)).copy(order='C')
        types = Trajectory.get_types_id_array(self.get_attribute('symbols'))
        params = [pos, vel, types, cel]
        atraj = pa.Trajectory(*params, pa.BoxFormat.CellVectors, wrapped, True)
        if save_reference and tskip==1:
            if wrapped:
                self.paTrajWrapped = atraj
            else:
                self.paTraj = atraj
        return atraj

    @staticmethod
    def get_type_mask_from_s(symbols):
        typesa = set(symbols)
        masks = {}
        types = []
        types_array = np.zeros(len(symbols), dtype=int)
        itype = 0
        for t in typesa:
            masks[t] = np.array(symbols) == t
            types.append(t)
            types_array[masks[t]] = itype
            itype += 1
        return types, types_array, masks

    def get_poscar(self, every=1, start=0, desc=''):
        """
        return array of poscar file strings.
        """
        poscars = []
        types, types_array, masks = Trajectory.get_type_mask_from_s(self.symbols)
        pos = self.get_array("positions")
        vel = self.get_array("velocities")
        cel = self.get_array("cells")
        nt = {}
        for it in types:
            nt[it] = np.sum(masks[it])
        for i in range(start, self.numsteps, every):
            c = cel[i]
            p = f'{desc}\n1.0\n'
            for idim in range(3):
                p += f'{c[idim, 0]} {c[idim, 1]} {c[idim, 2]}\n'
            for it in types:
                p += it + ' '
            p += '\n'
            for it in types:
                p += f'{nt[it]} '
            p += '\nCartesian\n'
            for it in types:
                ps = pos[i, masks[it]]
                for iatom in range(nt[it]):
                    p += f'{ps[iatom, 0]} {ps[iatom, 1]} {ps[iatom, 2]} \n'
            p += 'Cartesian\n'
            for it in types:
                v = vel[i, masks[it]]
                for iatom in range(nt[it]):
                    p += f'{v[iatom, 0]} {v[iatom, 1]} {v[iatom, 2]} \n'

            poscars.append(p)
        return poscars

    def rotate_qr(self,wrap=False,save_reference=DEFAULT_SAVE_REFERENCE):
        """
        transform the cell time series of this object to a triangular matrix.
        Rotate position, velocities, forces and stresses accordingly
        """
        wrapped=self.get_analisi_traj(wrapped=wrap,save_reference=save_reference)
        q = wrapped.get_rotation_matrix()
        new_cell = Trajectory.analisi_cell2box(wrapped.get_box_copy())
        self.data['cells']=new_cell
        self.data['positions']=wrapped.get_positions_copy()
        if 'forces' in self.get_arraynames():
            rotated_forces = np.einsum('taj,tij -> tai', self.get_array('forces'), q)
            self.data['forces'] = rotated_forces
        if 'stresses' in self.get_arraynames():
            rotated_stress = np.einsum('tlm,tli,tmj -> tij', self.get_array('stresses'), q, q)
            self.data['stresses']=rotated_stress
        if 'velocities' in self.get_arraynames():
            self.data['velocities'] = wrapped.get_velocities_copy()

    def get_aiida_traj(self):
        import aiida
        res = aiida.orm.nodes.data.array.trajectory.TrajectoryData()
        res.set_attribute('symbols', ft.symbols)
        for name in self.KEYS :
            if name in self.get_arraynames():
                res.set_array(name, ft.get_array(name))
        for name in self.KEYS_OPT :
            if name in self.get_arraynames():
                res.set_array(name, ft.get_array(name))
        for name in (self.TIME_K,self.STEP_K):
            if name in self.get_arraynames():
                res.set_array(name, ft.get_array(name))
        return res

    def get_lammps_data(self, timestep=0, comment=None, symb_id=None):
        def write_cell(a):
            res = f'{a[0]} {a[0] + 2 * a[3]} xlo xhi\n{a[1]} {a[1] + 2 * a[4]} ylo yhi\n{a[2]} {a[2] + 2 * a[5]} zlo zhi\n'
            if len(a) > 6:
                res += f'{a[6]} {a[7]} {a[8]} xy xz yz\n'
            return res

        wrapped = self.get_analisi_traj(wrapped=False,save_reference=self.DEFAULT_SAVE_REFERENCE)
        cel = wrapped.get_box_copy()
        pos = wrapped.get_positions_copy()
        vel = wrapped.get_velocities_copy()

        if timestep >= pos.shape[0] or timestep < -pos.shape[0]:
            raise IndexError('out of range index')
        if comment is None:
            comment = f'timestep {timestep}'
        if symb_id is None:
            all_sym = list(set(self.symbols))
            symb_id = {s: str(i + 1) for i, s in enumerate(all_sym)}
        symbols = [symb_id[s] for s in self.symbols]
        nt = len(symb_id)
        nat = pos.shape[1]

        res = f'{comment}\n{nat} atoms\n{nt} atom types\n{write_cell(cel[timestep])}\nAtoms\n\n'
        for iatom, itype, xyz in zip(range(nat), symbols, pos[timestep]):
            res += f'{iatom + 1} {itype} {xyz[0]} {xyz[1]} {xyz[2]}\n'
        res += '\nVelocities\n\n'
        for iatom, xyz in zip(range(nat), vel[timestep]):
            res += f'{iatom + 1} {xyz[0]} {xyz[1]} {xyz[2]}\n'

        return res

    def write_lammps_binary(self,fname, start=0, end=-1):
        if end > self.numsteps:
            raise IndexError(f'end must be <= {self.numsteps}')
        if start < 0 or start > self.numsteps:
            raise IndexError(f'start must be in [0, self.numsteps]')
        if 0 <= end < start:
            raise IndexError('start > end !')
        wrapped = self.get_analisi_traj(wrapped=False,save_reference=self.DEFAULT_SAVE_REFERENCE)
        wrapped.write_lammps_binary(fname, start, end)