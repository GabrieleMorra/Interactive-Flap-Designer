import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib.patches as patches
from scipy.interpolate import interp1d
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import json
import csv

class AirfoilFlapDesigner:
    def __init__(self, root):
        self.root = root
        self.root.title("Interactive Airfoil Flap Designer")
        self.root.state('zoomed')
        
        # Initialize data
        self.airfoil_data = None
        self.x_upper = None
        self.y_upper = None
        self.x_lower = None
        self.y_lower = None
        
        # Design parameters
        self.params = {
            'x1': 0.6, 'x2': 0.675, 'Delta_X1': 0.05, 'x3': 0.685, 'FractionY': 0.15,
            'DeltaX4': 0.02, 'Delta_angolo': 10, 'x6': 0.76, 'Delta_Lip': 0.045,
            'DeltaTE_main': 0.004, 'TwistMain': 0.0, 'x7': 0.84, 'DeltaX9': 0.04,
            'DeltaY9': 0.95, 'DeltaX10': 0.006, 'DeltaY10': -0.25, 'DeltaX11': 0.05,
            'x12': 0.73, 'xh': 0.75, 'yh': -0.1, 'FlapDeflection': 0.0
        }
        
        # Control points and interactive elements
        self.slot_control_points = None
        self.flap_control_points = None
        self.hinge_point = None
        self.slot_artists = []
        self.flap_artists = []
        self.hinge_artist = None
        self.point_annotations = []
        self.dragging_point = None
        self.drag_offset = None

        # Panning 
        self.is_panning = False
        self.pan_start = None
        self.pan_xlim = None
        self.pan_ylim = None
        
        # Magnetic snap
        self.magnetic_snap = False  
        self.snap_tolerance = 0.0025  # Tolleranza per lo snap

        # Gap and Overlap calculations
        self.gap_value = 0.0
        self.overlap_value = 0.0
        self.gap_overlap_vars = {}

        # Results
        self.main_component = None
        self.flap_component = None
        self.slot_curve = None
        self.flap_curve = None
        
        # Plot elements
        self.slot_curve_line = None
        self.flap_curve_line = None
        self.main_line = None
        self.flap_line = None
        
        # # Preview window
        # self.preview_window = None
        
        self.setup_gui()
        
    def setup_gui(self):
        # Create main frame with three sections
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel for parameters
        left_frame = ttk.Frame(main_frame, width=300)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        left_frame.pack_propagate(False)
        
        # Right panel for plot
        right_frame = ttk.Frame(main_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        self.setup_parameter_panel(left_frame)
        self.setup_plot_panel(right_frame)


    def setup_parameter_panel(self, parent):
        # Title
        ttk.Label(parent, text="Parameters", font=('Arial', 12, 'bold')).pack(pady=(0, 10))
        
        # Scrollable frame
        canvas = tk.Canvas(parent, highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Parameter controls
        self.param_vars = {}
        
        # Group parameters
        groups = {
            "Main Shape": ['x1', 'x2', 'Delta_X1', 'x3', 'FractionY', 'DeltaX4', 'Delta_angolo', 'x6', 'Delta_Lip', 'DeltaTE_main'],
            "Flap Shape": ['x7', 'DeltaX9', 'DeltaY9', 'DeltaX10', 'DeltaY10', 'DeltaX11', 'x12'],
            "Hinge & Deflection": ['xh', 'yh', 'FlapDeflection', 'TwistMain'],
            "Flap Design": ['Gap', 'Overlap']
        }
        
        for group_name, param_list in groups.items():
            # Group label
            group_frame = ttk.LabelFrame(scrollable_frame, text=group_name, padding=5)
            group_frame.pack(fill=tk.X, pady=5)
            
            for param in param_list:
                if param in self.params:
                    frame = ttk.Frame(group_frame)
                    frame.pack(fill=tk.X, pady=2)
                    
                    ttk.Label(frame, text=f"{param}:", width=12).pack(side=tk.LEFT)
                    
                    var = tk.DoubleVar(value=self.params[param])
                    self.param_vars[param] = var
                    
                    entry = ttk.Entry(frame, textvariable=var, width=10)
                    entry.pack(side=tk.LEFT, padx=(5, 0))
                    entry.bind('<Return>', lambda e, p=param: self.update_parameter(p))
                    entry.bind('<FocusOut>', lambda e, p=param: self.update_parameter(p))

        # Gap and Overlap display section (READ-ONLY)
        gap_overlap_frame = ttk.LabelFrame(scrollable_frame, text="Gap & Overlap (Read-only)", padding=5)
        gap_overlap_frame.pack(fill=tk.X, pady=5)
        
        # Gap display
        gap_frame = ttk.Frame(gap_overlap_frame)
        gap_frame.pack(fill=tk.X, pady=2)
        ttk.Label(gap_frame, text="Gap:", width=12).pack(side=tk.LEFT)
        gap_var = tk.StringVar(value="0.000000")
        self.gap_overlap_vars['gap'] = gap_var
        gap_entry = ttk.Entry(gap_frame, textvariable=gap_var, width=12, state='readonly')
        gap_entry.pack(side=tk.LEFT, padx=(5, 0))
        
        # Overlap display
        overlap_frame = ttk.Frame(gap_overlap_frame)
        overlap_frame.pack(fill=tk.X, pady=2)
        ttk.Label(overlap_frame, text="Overlap:", width=12).pack(side=tk.LEFT)
        overlap_var = tk.StringVar(value="0.000000")
        self.gap_overlap_vars['overlap'] = overlap_var
        overlap_entry = ttk.Entry(overlap_frame, textvariable=overlap_var, width=12, state='readonly')
        overlap_entry.pack(side=tk.LEFT, padx=(5, 0))
        
        # Control buttons
        btn_frame = ttk.Frame(scrollable_frame)
        btn_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(btn_frame, text="Update All", command=self.update_all_parameters).pack(fill=tk.X, pady=2)
        # ttk.Button(btn_frame, text="Show Preview", command=self.show_flap_preview).pack(fill=tk.X, pady=2)
        # ttk.Button(btn_frame, text="Export CSV", command=self.export_csv_dialog).pack(fill=tk.X, pady=2)
        
    def setup_plot_panel(self, parent):
        # Top control panel
        control_frame = ttk.Frame(parent)
        control_frame.pack(fill=tk.X, pady=(0, 5))
        
        ttk.Button(control_frame, text="Load Airfoil", command=self.load_airfoil).pack(side=tk.LEFT, padx=2)
        ttk.Button(control_frame, text="Reset View", command=self.reset_view).pack(side=tk.LEFT, padx=2)
        self.magnet_button = ttk.Button(control_frame, text="🧲 Snap", command=self.toggle_magnetic_snap)
        self.magnet_button.pack(side=tk.LEFT, padx=2)
        self.magnet_style = ttk.Style()

        ttk.Button(control_frame, text="Save Config", command=self.save_parameters).pack(side=tk.LEFT, padx=2)
        ttk.Button(control_frame, text="Load Config", command=self.load_parameters).pack(side=tk.LEFT, padx=2)

        ttk.Button(control_frame, text="Export", command=self.export_dialog).pack(side=tk.LEFT, padx=2)

        
        self.status_label = ttk.Label(control_frame, text="Drag control points to modify shape")
        self.status_label.pack(side=tk.RIGHT, padx=10)
        
        # Plot frame
        plot_frame = ttk.Frame(parent)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        
        self.fig = Figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Bind mouse events
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect("scroll_event", self.on_mouse_scroll)

    def toggle_magnetic_snap(self):
        """Toggle magnetic snap on/off"""
        self.magnetic_snap = not self.magnetic_snap
        if self.magnetic_snap:
            # Creare uno style personalizzato per il bottone attivo
            self.magnet_style.configure("Active.TButton", background="lightgreen")
            self.magnet_button.config(style="Active.TButton")
            self.status_label.config(text="Magnetic snap ON - Points will snap to original airfoil")
        else:
            self.magnet_button.config(style="TButton")  # Stile default
            self.status_label.config(text="Magnetic snap OFF")

    def find_nearest_airfoil_point(self, x, y):
        """Find nearest point on original airfoil curve (interpolated)"""
        if self.airfoil_data is None:
            return x, y
        
        # Creare interpolazione ad alta risoluzione della curva originale
        x_curve = self.airfoil_data[:, 0]
        y_curve = self.airfoil_data[:, 1]
        
        # Creare punti interpolati ad alta densità
        t_original = np.linspace(0, 1, len(x_curve))
        t_interp = np.linspace(0, 1, len(x_curve) * 20)  # 20x più punti
        
        from scipy.interpolate import interp1d
        f_x = interp1d(t_original, x_curve, kind='linear')
        f_y = interp1d(t_original, y_curve, kind='linear')
        
        x_interp = f_x(t_interp)
        y_interp = f_y(t_interp)
        
        # Calcola distanze da tutti i punti interpolati
        distances = np.sqrt((x_interp - x)**2 + (y_interp - y)**2)
        min_dist_idx = np.argmin(distances)
        min_distance = distances[min_dist_idx]
        
        # Se la distanza minima è entro la tolleranza, snap al punto interpolato
        if min_distance < self.snap_tolerance:
            return x_interp[min_dist_idx], y_interp[min_dist_idx]
        else:
            return x, y
        
    def export_dialog(self):
        """Open export dialog window"""
        if self.main_component is None:
            messagebox.showwarning("Warning", "No data to export. Please generate a flap first!")
            return
        
        export_window = tk.Toplevel(self.root)
        export_window.title("Export Settings")
        export_window.geometry("300x250")
        export_window.resizable(False, False)
        
        # Variables
        self.export_main = tk.BooleanVar(value=True)
        self.export_flap = tk.BooleanVar(value=True)
        self.main_points = tk.IntVar(value=150)
        self.flap_points = tk.IntVar(value=75)
        
        # Main frame
        frame = ttk.Frame(export_window, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Export options
        ttk.Label(frame, text="Export Options:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(0,5))
        
        main_check = ttk.Checkbutton(frame, text="Export Main Component", variable=self.export_main,
                                    command=self.update_export_options)
        main_check.pack(anchor='w')
        
        flap_check = ttk.Checkbutton(frame, text="Export Flap Component", variable=self.export_flap,
                                    command=self.update_export_options)
        flap_check.pack(anchor='w', pady=(0,10))
        
        # Points selection
        ttk.Label(frame, text="Number of Points:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(10,5))
        
        # Main points
        main_frame = ttk.Frame(frame)
        main_frame.pack(fill='x', pady=2)
        self.main_label = ttk.Label(main_frame, text="Main points:")
        self.main_label.pack(side='left')
        self.main_entry = ttk.Entry(main_frame, textvariable=self.main_points, width=8)
        self.main_entry.pack(side='right')
        
        # Flap points
        flap_frame = ttk.Frame(frame)
        flap_frame.pack(fill='x', pady=2)
        self.flap_label = ttk.Label(flap_frame, text="Flap points:")
        self.flap_label.pack(side='left')
        self.flap_entry = ttk.Entry(flap_frame, textvariable=self.flap_points, width=8)
        self.flap_entry.pack(side='right')
        
        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', pady=(20,0))
        
        ttk.Button(btn_frame, text="Export", command=lambda: self.perform_export(export_window)).pack(side='right', padx=(5,0))
        ttk.Button(btn_frame, text="Cancel", command=export_window.destroy).pack(side='right')
        
        self.update_export_options()

    def calculate_gap_overlap(self):
        """Calculate gap and overlap between main component and flap"""
        if self.main_component is None or self.flap_component is None:
            self.gap_value = 0.0
            self.overlap_value = 0.0
            return
        
        try:
            # Get deflected flap for calculations
            deflected_flap = self.apply_flap_deflection()
            
            # Find main component trailing edge (maximum x coordinate)
            main_te_idx = np.argmax(self.main_component[:, 0])
            main_te_point = self.main_component[main_te_idx]
            
            # Find flap leading edge (minimum x coordinate)
            flap_le_idx = np.argmin(deflected_flap[:, 0])
            flap_le_point = deflected_flap[flap_le_idx]
            
            # === OVERLAP CALCULATION ===
            # Overlap: distanza orizzontale tra punto più dietro del main e punto più avanti del flap
            # Positivo se il flap sta a x più basse del main
            x_main_max = np.max(self.main_component[:, 0])  # Punto più dietro del main
            x_flap_min = np.min(deflected_flap[:, 0])       # Punto più avanti del flap
            
            self.overlap_value = x_main_max - x_flap_min
            
            # === GAP CALCULATION ===
            # Gap: distanza minima tra trailing edge del main e qualsiasi punto del flap
            min_distance = float('inf')
            
            # Calcola la distanza minima tra il trailing edge e tutti i punti del flap
            for flap_point in deflected_flap:
                distance = np.sqrt((main_te_point[0] - flap_point[0])**2 + 
                                (main_te_point[1] - flap_point[1])**2)
                if distance < min_distance:
                    min_distance = distance
            
            self.gap_value = min_distance
            
            # Se c'è overlap (flap avanza oltre il main), il gap potrebbe essere zero
            # o molto piccolo se le superfici si sovrappongono
            if self.overlap_value > 0:
                # Verifica se c'è effettiva sovrapposizione geometrica
                # che potrebbe rendere il gap praticamente zero
                pass  # Il gap rimane quello calcolato come distanza minima
                
        except Exception as e:
            print(f"Error calculating gap/overlap: {str(e)}")
            self.gap_value = 0.0
            self.overlap_value = 0.0
        
        # Update GUI variables
        if hasattr(self, 'gap_overlap_vars'):
            if 'gap' in self.gap_overlap_vars:
                self.gap_overlap_vars['gap'].set(f"{self.gap_value:.6f}")
            if 'overlap' in self.gap_overlap_vars:
                self.gap_overlap_vars['overlap'].set(f"{self.overlap_value:.6f}")

    def update_export_options(self):
        """Update export options based on checkboxes"""
        main_enabled = self.export_main.get()
        flap_enabled = self.export_flap.get()
        
        # Enable/disable entries
        self.main_entry.config(state='normal' if main_enabled else 'disabled')
        self.flap_entry.config(state='normal' if flap_enabled else 'disabled')
        self.main_label.config(foreground='black' if main_enabled else 'gray')
        self.flap_label.config(foreground='black' if flap_enabled else 'gray')

    def perform_export(self, export_window):
        """Perform the actual export"""
        if not self.export_main.get() and not self.export_flap.get():
            messagebox.showwarning("Warning", "Please select at least one component to export!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export to CSV",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                self.export_to_csv(filename)
                export_window.destroy()
                self.status_label.config(text="Export completed successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")

    def export_to_csv(self, filename):
        """Export components to CSV file"""
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            
            if self.export_main.get():
                writer.writerow("")
                writer.writerow(['Main coords (x,y,z)'])
                main_data = self.resample_component(self.main_component, self.main_points.get())
                for point in main_data:
                    writer.writerow([point[0], point[1], 0.0])
            
            if self.export_flap.get():
                writer.writerow("")
                writer.writerow(['Flap coords (x,y,z)'])
                flap_deflected = self.apply_flap_deflection()
                flap_data = self.resample_component(flap_deflected, self.flap_points.get())
                for point in flap_data:
                    writer.writerow([point[0], point[1], 0.0])

    def apply_flap_deflection(self):
        """Apply rigid rotation to flap around hinge point"""
        if self.flap_component is None or abs(self.params['FlapDeflection']) < 1e-6:
            return self.flap_component.copy() if self.flap_component is not None else None
        
        # Get hinge point and deflection angle
        xh, yh = self.params['xh'], self.params['yh']
        delta = np.deg2rad(self.params['FlapDeflection'])  # Positive = nose down
        
        # Translate flap so hinge is at origin
        x_trans = self.flap_component[:, 0] - xh
        y_trans = self.flap_component[:, 1] - yh
        
        # Apply rotation matrix (clockwise for positive deflection)
        cos_delta = np.cos(-delta)  # Negative for clockwise rotation
        sin_delta = np.sin(-delta)
        
        x_rot = x_trans * cos_delta - y_trans * sin_delta
        y_rot = x_trans * sin_delta + y_trans * cos_delta
        
        # Translate back
        deflected_flap = np.column_stack([x_rot + xh, y_rot + yh])
        
        return deflected_flap

    def resample_component(self, component, n_points):
        """Resample component to specified number of points, starting from trailing edge"""
        if component is None or len(component) == 0:
            return np.array([])
        
        # Find trailing edge (maximum x coordinate)
        te_idx = np.argmax(component[:, 0])
        
        # Reorder starting from trailing edge: TE -> upper surface -> LE -> lower surface -> TE
        reordered = np.concatenate([component[te_idx:], component[:te_idx+1]])
        
        # Create parameter t for interpolation
        distances = np.concatenate([[0], np.cumsum(np.sqrt(np.diff(reordered[:, 0])**2 + np.diff(reordered[:, 1])**2))])
        t_original = distances / distances[-1]
        
        # New parameter values for resampling
        t_new = np.linspace(0, 1, n_points)
        
        # Interpolate
        from scipy.interpolate import interp1d
        fx = interp1d(t_original, reordered[:, 0], kind='linear')
        fy = interp1d(t_original, reordered[:, 1], kind='linear')
        
        resampled = np.column_stack([fx(t_new), fy(t_new)])
        return resampled

    
    def update_parameter(self, param_name):
        """Update a single parameter and regenerate"""
        try:
            new_value = self.param_vars[param_name].get()
            self.params[param_name] = new_value
            
            # Update corresponding control point if it exists
            if param_name == 'xh' and self.hinge_point is not None:
                self.hinge_point[0] = new_value
            elif param_name == 'yh' and self.hinge_point is not None:
                self.hinge_point[1] = new_value
            
            self.generate_flap()
        except Exception as e:
            self.status_label.config(text=f"Error updating {param_name}: {str(e)}")
    
    def update_all_parameters(self):
        """Update all parameters and regenerate"""
        try:
            # Salva i limiti correnti
            current_xlim = self.ax.get_xlim()
            current_ylim = self.ax.get_ylim()
            
            for param, var in self.param_vars.items():
                self.params[param] = var.get()
            self.generate_flap()
            
            # Ripristina i limiti
            self.ax.set_xlim(current_xlim)
            self.ax.set_ylim(current_ylim)
            self.canvas.draw()
        except Exception as e:
            self.status_label.config(text=f"Error updating parameters: {str(e)}")

    def update_parameters_from_sketch(self):
        """Update all parameters based on control points moved in the sketch"""
        try:
            if self.slot_control_points is None or self.flap_control_points is None:
                return
            
            p = self.params
            f_upper = interp1d(self.x_upper, self.y_upper, kind='linear', bounds_error=False, fill_value='extrapolate')
            f_lower = interp1d(self.x_lower, self.y_lower, kind='linear', bounds_error=False, fill_value='extrapolate')
            
            # Extract control points
            x1, y1 = self.slot_control_points[0]
            x2, y2 = self.slot_control_points[1] 
            x3, y3 = self.slot_control_points[2]
            x4, y4 = self.slot_control_points[3]
            x5, y5 = self.slot_control_points[4]
            x6, y6 = self.slot_control_points[5]
            
            x7, y7 = self.flap_control_points[0]
            x8, y8 = self.flap_control_points[1]
            x9, y9 = self.flap_control_points[2]
            x10, y10 = self.flap_control_points[3]
            x11, y11 = self.flap_control_points[4]
            x12, y12 = self.flap_control_points[5]
            
            # Direct coordinate updates
            p['x1'] = x1
            p['x2'] = x2
            p['x3'] = x3
            p['x6'] = x6
            p['x7'] = x7
            p['x12'] = x12
            p['xh'] = self.hinge_point[0]
            p['yh'] = self.hinge_point[1]
            
            # Delta_X1: calcolo dalla geometria attuale
            y1_orig = f_lower(x1)
            # Calcola la pendenza della curva originale in x1
            dx = 0.001
            y1_left = f_lower(x1 - dx)
            slope_orig = (y1_orig - y1_left) / dx
            
            # Calcola la pendenza attuale dal punto 1 al punto 2
            if x2 != x1:
                slope_actual = (y2 - y1) / (x2 - x1)
                # Delta_X1 è la distanza necessaria per ottenere questa pendenza
                if slope_orig != 0:
                    p['Delta_X1'] = abs((slope_actual - slope_orig) / slope_orig * 0.05)
                else:
                    p['Delta_X1'] = abs(x1 - x2) * 0.5
            
            # FractionY: frazione tra slot line e upper surface
            y3_upper = f_upper(x3)
            if abs(y3_upper - y2) > 1e-6:
                p['FractionY'] = max(0, min(1, (y3 - y2) / (y3_upper - y2)))
            
            # DeltaX4
            p['DeltaX4'] = abs(x3 - x4)
            
            # Delta_angolo: CORREZIONE QUI
            if x3 != x2 and x4 != x3:
                # Calcola angolo_3 come nel codice originale (generate_slot_control_points)
                angolo_3 = 90 - 57.3 * np.arctan((y3 - y2) / (x3 - x2))
                
                # Calcola angolo_tot dal vettore 3->4
                slope_34 = (y4 - y3) / (x4 - x3)
                angolo_tot = np.degrees(np.arctan(slope_34))
                
                # Delta_angolo = angolo_3 - angolo_tot (come nella definizione originale)
                p['Delta_angolo'] = angolo_3 - angolo_tot
            
            # Delta_Lip
            p['Delta_Lip'] = abs(x6 - x5)
            
            # DeltaTE_main
            y6_upper = f_upper(x6)
            p['DeltaTE_main'] = abs(y6_upper - y6)
            
            # Flap parameters
            p['DeltaX9'] = x9 - x3
            
            # DeltaY9: normalizzato rispetto alla distanza verticale
            vdist = abs(f_upper(x6) - f_lower(x6))
            if vdist > 1e-6:
                p['DeltaY9'] = (y9 - y3) / vdist
            
            # DeltaX10
            p['DeltaX10'] = x10 - x9
            
            # DeltaY10: normalizzato
            if vdist > 1e-6:
                p['DeltaY10'] = (y10 - y9) / vdist
            
            # DeltaX11
            p['DeltaX11'] = abs(x12 - x11)
            
            # Update parameter GUI fields
            for param, var in self.param_vars.items():
                if param in p:
                    var.set(round(p[param], 6))
                    
        except Exception as e:
            self.status_label.config(text=f"Error updating parameters: {str(e)}")

    def load_airfoil(self):
        """Load airfoil data from file"""
        filename = filedialog.askopenfilename(
            title="Select Airfoil File",
            filetypes=[("Text files", "*.txt"), ("DAT files", "*.dat"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                self.airfoil_data = np.loadtxt(filename)
                self.split_airfoil()
                self.plot_airfoil()
                self.generate_flap()
                self.status_label.config(text="Airfoil loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Error loading airfoil: {str(e)}")
    
    def split_airfoil(self):
        """Split airfoil into upper and lower surfaces"""
        if self.airfoil_data is None:
            return
        x_le_idx = np.argmin(self.airfoil_data[:, 0])
        self.x_upper = self.airfoil_data[:x_le_idx+1, 0]
        self.y_upper = self.airfoil_data[:x_le_idx+1, 1]
        self.x_lower = self.airfoil_data[x_le_idx:, 0]
        self.y_lower = self.airfoil_data[x_le_idx:, 1]
    
    def bezier_curve(self, control_points, n_points=60):
        """Generate Bezier curve from control points"""
        n = len(control_points) - 1
        t = np.linspace(0, 1, n_points)
        curve = np.zeros((n_points, 2))
        
        for i in range(n_points):
            for j in range(n + 1):
                binom_coeff = math.comb(n, j)
                bernstein = binom_coeff * (t[i] ** j) * ((1 - t[i]) ** (n - j))
                curve[i] += bernstein * control_points[j]
        
        return curve
    
    def generate_flap(self):
        """Generate flap based on current parameters"""
        if self.airfoil_data is None:
            self.status_label.config(text="Please load an airfoil first!")
            return
        
        try:
            self.generate_slot_control_points()
            self.generate_flap_control_points()
            self.hinge_point = np.array([self.params['xh'], self.params['yh']])
            self.slot_curve = self.bezier_curve(self.slot_control_points)
            self.flap_curve = self.bezier_curve(self.flap_control_points)
            self.assemble_components()
            self.plot_results()
            self.status_label.config(text="Flap generated! Drag control points to modify")
            self.calculate_gap_overlap()
        except Exception as e:
            self.status_label.config(text=f"Error generating flap: {str(e)}")
        
    
    def generate_slot_control_points(self):
        """Generate slot control points"""
        p = self.params
        f_upper = interp1d(self.x_upper, self.y_upper, kind='linear', bounds_error=False, fill_value='extrapolate')
        f_lower = interp1d(self.x_lower, self.y_lower, kind='linear', bounds_error=False, fill_value='extrapolate')
        
        # Point calculations (same as original)
        x1, y1 = p['x1'], f_lower(p['x1'])
        x2p = x1 - p['Delta_X1']
        y2p = f_lower(x2p)
        slope_lower = np.arctan((y2p - y1) / (x2p - x1))
        x2 = p['x2']
        y2 = y1 + slope_lower * (x2 - x1)
        
        x3 = p['x3']
        y3p = f_upper(x3)
        delta_y3 = p['FractionY'] * (y3p - y2)
        y3 = y2 + delta_y3
        
        x4 = x3 - p['DeltaX4']
        angolo_3 = 90 - 57.3 * np.arctan((y3 - y2) / (x3 - x2))
        angolo_tot = angolo_3 - p['Delta_angolo']
        slope_4 = np.tan(np.deg2rad(angolo_tot))
        y4 = y3 + slope_4 * (x4 - x3)
        
        x6 = p['x6']
        y6_upper = f_upper(x6)
        y6 = y6_upper - p['DeltaTE_main']
        
        x6p = x6 + 0.05
        y6p = f_upper(x6p)
        slope_p6 = -np.arctan((y6p - y6) / (x6p - x6))
        x5 = x6 - p['Delta_Lip']
        y5 = y6 + slope_p6 * (x6 - x5)
        
        self.slot_control_points = np.array([[x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5], [x6, y6]])
    
    def generate_flap_control_points(self):
        """Generate flap control points"""
        p = self.params
        f_upper = interp1d(self.x_upper, self.y_upper, kind='linear', bounds_error=False, fill_value='extrapolate')
        f_lower = interp1d(self.x_lower, self.y_lower, kind='linear', bounds_error=False, fill_value='extrapolate')
        
        x7, y7 = p['x7'], f_upper(p['x7'])
        x8, y8 = p['x6'], f_upper(p['x6'])
        
        x9 = self.slot_control_points[2, 0] + p['DeltaX9']
        vdist_9 = abs(f_upper(p['x6']) - f_lower(p['x6']))
        y9 = self.slot_control_points[2, 1] + p['DeltaY9'] * vdist_9
        
        x10 = x9 + p['DeltaX10']
        y10 = y9 + p['DeltaY10'] * vdist_9
        
        x12, y12 = p['x12'], f_lower(p['x12'])
        x11 = x12 - p['DeltaX11']
        x12_bis = x12 + 0.00001
        y12_bis = f_lower(x12_bis)
        slope_y11 = np.arctan((y12_bis - y12) / (x12_bis - x12))
        y11 = y12 - slope_y11 * (x12 - x11)
        
        self.flap_control_points = np.array([[x7, y7], [x8, y8], [x9, y9], [x10, y10], [x11, y11], [x12, y12]])
    
    def assemble_components(self):
        """Assemble main and flap components"""
        p = self.params
        
        ix6 = np.argmin(np.abs(p['x6'] - self.x_upper))
        ix1 = np.argmin(np.abs(p['x1'] - self.x_lower))
        ix7 = np.argmin(np.abs(p['x7'] - self.x_upper))
        ix12 = np.argmin(np.abs(p['x12'] - self.x_lower))
        
        # Main component
        main_x = np.concatenate([self.x_upper[ix6:], self.x_lower[:ix1], self.slot_curve[:, 0]])
        main_y = np.concatenate([self.y_upper[ix6:], self.y_lower[:ix1], self.slot_curve[:, 1]])
        
        cos_twist = np.cos(np.deg2rad(p['TwistMain']))
        sin_twist = np.sin(np.deg2rad(p['TwistMain']))
        self.main_component = np.column_stack([
            main_x * cos_twist - main_y * sin_twist,
            main_x * sin_twist + main_y * cos_twist
        ])
        
        # Flap component (undeflected)
        flap_x = np.concatenate([self.x_upper[:ix7], self.flap_curve[:, 0], self.x_lower[ix12+1:]])
        flap_y = np.concatenate([self.y_upper[:ix7], self.flap_curve[:, 1], self.y_lower[ix12+1:]])
        
        # Apply main rotation only for display
        self.flap_component = np.column_stack([
            flap_x * cos_twist - flap_y * sin_twist,
            flap_x * sin_twist + flap_y * cos_twist
        ])
    
    def get_deflected_flap(self):
        """Get flap component with deflection applied"""
        p = self.params
        
        # Apply deflection
        xh, yh = p['xh'], p['yh']
        delta = np.deg2rad(p['FlapDeflection'])
        
        x_trans = self.flap_component[:, 0] - xh
        y_trans = self.flap_component[:, 1] - yh
        
        cos_delta = np.cos(-delta)
        sin_delta = np.sin(-delta)
        x_rot = x_trans * cos_delta - y_trans * sin_delta
        y_rot = x_trans * sin_delta + y_trans * cos_delta
        
        return np.column_stack([x_rot + xh, y_rot + yh])
    
    def plot_airfoil(self):
        """Plot the original airfoil"""
        self.ax.clear()
        if self.airfoil_data is not None:
            self.ax.plot(self.airfoil_data[:, 0], self.airfoil_data[:, 1], 'k-', linewidth=1, alpha=0.5, label='Original')
            self.ax.set_xlabel('x/c')
            self.ax.set_ylabel('y/c')
            self.ax.grid(True, alpha=0.3)
            self.ax.axis('equal')
            self.ax.legend()
        self.canvas.draw()
    
    def plot_results(self):
        """Plot results without deflected flap"""
        self.ax.clear()
        
        if self.airfoil_data is not None:
            self.ax.plot(self.airfoil_data[:, 0], self.airfoil_data[:, 1], 'k-', linewidth=1, alpha=0.3, label='Original')
        
        if self.main_component is not None:
            self.main_line = self.ax.plot(self.main_component[:, 0], self.main_component[:, 1], 'b-', linewidth=2.5, label='Main')[0]
        
        if self.flap_component is not None:
            # Plot undeflected flap
            self.flap_line = self.ax.plot(self.flap_component[:, 0], self.flap_component[:, 1], 'r--', linewidth=1.5, alpha=0.5, label='Flap (undeflected)')[0]
            
            # Plot deflected flap if deflection is not zero
            if abs(self.params['FlapDeflection']) > 1e-6:
                deflected_flap = self.apply_flap_deflection()
                self.deflected_flap_line = self.ax.plot(deflected_flap[:, 0], deflected_flap[:, 1], 'r-', linewidth=2.5, label=f'Flap (deflected {self.params["flap_deflection"]:.1f}°)')[0]

        if self.slot_curve is not None:
            self.slot_curve_line = self.ax.plot(self.slot_curve[:, 0], self.slot_curve[:, 1], 'b:', linewidth=1.5, alpha=0.7, label='Main Curve')[0]
        
        if self.flap_curve is not None:
            self.flap_curve_line = self.ax.plot(self.flap_curve[:, 0], self.flap_curve[:, 1], 'r:', linewidth=1.5, alpha=0.7, label='Flap Curve')[0]

        self.create_interactive_points()
        
        self.ax.set_xlabel('x/c')
        self.ax.set_ylabel('y/c')
        self.ax.grid(True, alpha=0.3)
        self.ax.axis('equal')
        self.ax.legend()
        self.ax.set_title('Interactive Airfoil Flap Designer')
        
        self.canvas.draw()
    
    def create_interactive_points(self):
        """Create interactive control points with better hit detection"""
        # Clear existing artists and annotations
        for artist in self.slot_artists + self.flap_artists:
            if artist in self.ax.collections:
                artist.remove()
        if self.hinge_artist and self.hinge_artist in self.ax.collections:
            self.hinge_artist.remove()
        
        # Clear existing control lines
        if hasattr(self, 'slot_control_lines'):
            for line in self.slot_control_lines:
                line.remove()
        if hasattr(self, 'flap_control_lines'):
            for line in self.flap_control_lines:
                line.remove()
        
        # RIMUOVERE le annotazioni esistenti
        for ann in self.point_annotations:
            ann.remove()
        
        self.slot_artists = []
        self.flap_artists = []
        self.slot_control_lines = []
        self.flap_control_lines = []
        self.point_annotations = []  # Mantenere vuoto
        
        # Create slot control lines
        if self.slot_control_points is not None:
            for i in range(len(self.slot_control_points) - 1):
                x_data = [self.slot_control_points[i][0], self.slot_control_points[i+1][0]]
                y_data = [self.slot_control_points[i][1], self.slot_control_points[i+1][1]]
                line = self.ax.plot(x_data, y_data, 'b--', linewidth=0.8, alpha=0.6, zorder=5)[0]
                self.slot_control_lines.append(line)
        
        # Create flap control lines
        if self.flap_control_points is not None:
            for i in range(len(self.flap_control_points) - 1):
                x_data = [self.flap_control_points[i][0], self.flap_control_points[i+1][0]]
                y_data = [self.flap_control_points[i][1], self.flap_control_points[i+1][1]]
                line = self.ax.plot(x_data, y_data, 'r--', linewidth=0.8, alpha=0.6, zorder=5)[0]
                self.flap_control_lines.append(line)
        
        # Create slot control points
        if self.slot_control_points is not None:
            for i, point in enumerate(self.slot_control_points):
                artist = self.ax.scatter(point[0], point[1], c='blue', s=80, marker='o', 
                                    edgecolor='darkblue', linewidth=2, zorder=15, alpha=0.9)
                self.slot_artists.append(artist)
        
        # Create flap control points
        if self.flap_control_points is not None:
            for i, point in enumerate(self.flap_control_points):
                artist = self.ax.scatter(point[0], point[1], c='red', s=80, marker='s', 
                                    edgecolor='darkred', linewidth=2, zorder=15, alpha=0.9)
                self.flap_artists.append(artist)
        
        # Create hinge point
        if self.hinge_point is not None:
            self.hinge_artist = self.ax.scatter(self.hinge_point[0], self.hinge_point[1], 
                                            c='green', s=100, marker='D', 
                                            edgecolor='darkgreen', linewidth=2, zorder=15, alpha=0.9)

    def update_curves_and_components(self):
        """Update curves and components"""
        try:
            self.slot_curve = self.bezier_curve(self.slot_control_points)
            self.flap_curve = self.bezier_curve(self.flap_control_points)
            
            # Update parameter values
            self.params['xh'] = self.hinge_point[0]
            self.params['yh'] = self.hinge_point[1]
            self.param_vars['xh'].set(self.hinge_point[0])
            self.param_vars['yh'].set(self.hinge_point[1])
            
            self.assemble_components()
            
            # Update plot lines
            if self.main_line:
                self.main_line.set_data(self.main_component[:, 0], self.main_component[:, 1])
            if self.flap_line:
                self.flap_line.set_data(self.flap_component[:, 0], self.flap_component[:, 1])
            if self.slot_curve_line:
                self.slot_curve_line.set_data(self.slot_curve[:, 0], self.slot_curve[:, 1])
            if self.flap_curve_line:
                self.flap_curve_line.set_data(self.flap_curve[:, 0], self.flap_curve[:, 1])
            
            # Update control lines
            if hasattr(self, 'slot_control_lines') and self.slot_control_points is not None:
                for i, line in enumerate(self.slot_control_lines):
                    if i < len(self.slot_control_points) - 1:
                        x_data = [self.slot_control_points[i][0], self.slot_control_points[i+1][0]]
                        y_data = [self.slot_control_points[i][1], self.slot_control_points[i+1][1]]
                        line.set_data(x_data, y_data)
            
            if hasattr(self, 'flap_control_lines') and self.flap_control_points is not None:
                for i, line in enumerate(self.flap_control_lines):
                    if i < len(self.flap_control_points) - 1:
                        x_data = [self.flap_control_points[i][0], self.flap_control_points[i+1][0]]
                        y_data = [self.flap_control_points[i][1], self.flap_control_points[i+1][1]]
                        line.set_data(x_data, y_data)
            
            # Update deflected flap line if it exists
            if hasattr(self, 'deflected_flap_line') and abs(self.params['FlapDeflection']) > 1e-6:
                deflected_flap = self.apply_flap_deflection()
                self.deflected_flap_line.set_data(deflected_flap[:, 0], deflected_flap[:, 1])
                # Update label
                self.deflected_flap_line.set_label(f'Flap (deflected {self.params["flap_deflection"]:.1f}°)')
            elif hasattr(self, 'deflected_flap_line') and abs(self.params['FlapDeflection']) <= 1e-6:
                # Remove deflected flap line if deflection is zero
                self.deflected_flap_line.remove()
                delattr(self, 'deflected_flap_line')
            elif not hasattr(self, 'deflected_flap_line') and abs(self.params['FlapDeflection']) > 1e-6:
                # Create deflected flap line if it doesn't exist
                deflected_flap = self.apply_flap_deflection()
                self.deflected_flap_line = self.ax.plot(deflected_flap[:, 0], deflected_flap[:, 1], 'r-', linewidth=2.5, 
                                                    label=f'Flap (deflected {self.params["flap_deflection"]:.1f}°)')[0]
            
            # Update point annotations
            for i, ann in enumerate(self.point_annotations):
                if i < len(self.slot_control_points):
                    ann.set_position(self.slot_control_points[i])
                elif i < len(self.slot_control_points) + len(self.flap_control_points):
                    idx = i - len(self.slot_control_points)
                    ann.set_position(self.flap_control_points[idx])
                else:
                    ann.set_position(self.hinge_point)

            self.calculate_gap_overlap()

            self.canvas.draw_idle()
            
        except Exception as e:
            self.status_label.config(text=f"Error updating: {str(e)}")

    def on_mouse_scroll(self, event):
        """Handle mouse scroll for zooming"""
        if event.inaxes != self.ax:
            return
        
        scale_factor = 1.1 if event.step > 0 else 0.9
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        
        x_center = event.xdata
        y_center = event.ydata
        
        new_xlim = [(x - x_center) * scale_factor + x_center for x in xlim]
        new_ylim = [(y - y_center) * scale_factor + y_center for y in ylim]
        
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        self.canvas.draw_idle()

    def on_mouse_press(self, event):
        """Handle mouse press with improved hit detection"""
        if event.inaxes != self.ax:
            return
        
        # if left click:
        if event.button == 1:
            tolerance = 0.008  # Ridotta da 0.015 per maggiore precisione
            
            # Raccogliere tutti i punti con le loro distanze
            candidates = []
            
            # Check slot control points
            if self.slot_control_points is not None:
                for i, point in enumerate(self.slot_control_points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance:
                        candidates.append((distance, 'slot', i, point))
            
            # Check flap control points
            if self.flap_control_points is not None:
                for i, point in enumerate(self.flap_control_points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance:
                        candidates.append((distance, 'flap', i, point))
            
            # Check hinge point
            if self.hinge_point is not None:
                distance = np.sqrt((event.xdata - self.hinge_point[0])**2 + (event.ydata - self.hinge_point[1])**2)
                if distance < tolerance:
                    candidates.append((distance, 'hinge', 0, self.hinge_point))
            
            # Selezionare il punto più vicino
            if candidates:
                # Ordinare per distanza (più vicino prima)
                candidates.sort(key=lambda x: x[0])
                distance, point_type, point_idx, point = candidates[0]
                
                self.dragging_point = (point_type, point_idx)
                self.drag_offset = (event.xdata - point[0], event.ydata - point[1])
                
                if point_type == 'slot':
                    self.status_label.config(text=f"Dragging Main Point {point_idx+1}")
                elif point_type == 'flap':
                    self.status_label.config(text=f"Dragging Flap Point {point_idx+1}")
                elif point_type == 'hinge':
                    self.status_label.config(text="Dragging Hinge Point")
                return
            
        elif event.button == 3:
            # Right click -> pan
            self.is_panning = True
            self.pan_start = (event.xdata, event.ydata)
            self.pan_xlim = self.ax.get_xlim()
            self.pan_ylim = self.ax.get_ylim()
            self.status_label.config(text="Panning view")

    def on_mouse_move(self, event):
        """Handle mouse move for point dragging"""
        if (self.dragging_point is None and self.pan_start is None) or event.inaxes != self.ax:
            return
        
        if self.dragging_point is not None and event.button == 1:
            curve_type, point_idx = self.dragging_point
            new_x = event.xdata - self.drag_offset[0]
            new_y = event.ydata - self.drag_offset[1]
            
            # APPLICARE SNAP MAGNETICO SE ATTIVO
            if self.magnetic_snap:
                new_x, new_y = self.find_nearest_airfoil_point(new_x, new_y)
            
            # Update control point
            if curve_type == 'slot':
                self.slot_control_points[point_idx] = [new_x, new_y]
                # Update corresponding scatter plot
                self.slot_artists[point_idx].set_offsets([[new_x, new_y]])
            elif curve_type == 'flap':
                self.flap_control_points[point_idx] = [new_x, new_y]
                # Update corresponding scatter plot
                self.flap_artists[point_idx].set_offsets([[new_x, new_y]])
            elif curve_type == 'hinge':
                self.hinge_point = np.array([new_x, new_y])
                # Update corresponding scatter plot
                self.hinge_artist.set_offsets([[new_x, new_y]])
            
            # Update curves and components
            self.update_curves_and_components()

        elif hasattr(self, 'pan_start') and self.pan_start is not None and self.is_panning and self.pan_xlim is not None and self.pan_ylim is not None:
            # MODIFICARE per pan fluido: calcolare offset in tempo reale
            dx = event.xdata - self.pan_start[0]
            dy = event.ydata - self.pan_start[1]

            new_xlim = [x - dx for x in self.pan_xlim]
            new_ylim = [y - dy for y in self.pan_ylim]

            self.ax.set_xlim(new_xlim)
            self.ax.set_ylim(new_ylim)


    def on_mouse_release(self, event):
        """Handle mouse release"""
        if self.dragging_point is not None:
            self.status_label.config(text="Drag control points to modify shape")
        self.dragging_point = None
        self.drag_offset = None

        if self.is_panning:
            self.canvas.draw_idle()
            self.is_panning = False
            self.pan_start = None  # AGGIUNGERE questa riga
            self.pan_xlim = None
            self.pan_ylim = None
            self.status_label.config(text="Drag control points to modify shape")
        
        self.update_parameters_from_sketch()
    
    def reset_view(self):
        """Reset the plot view"""
        self.ax.relim()
        self.ax.autoscale()
        self.canvas.draw()
    
    def save_parameters(self):
        """Save current parameters and control points to file"""
        filename = filedialog.asksaveasfilename(
            title="Save Configuration",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                config = {
                    'parameters': self.params.copy(),
                    'slot_control_points': self.slot_control_points.tolist() if self.slot_control_points is not None else None,
                    'flap_control_points': self.flap_control_points.tolist() if self.flap_control_points is not None else None,
                    'hinge_point': self.hinge_point.tolist() if self.hinge_point is not None else None,
                    'magnetic_snap': self.magnetic_snap
                }
                
                with open(filename, 'w') as f:
                    json.dump(config, f, indent=2)
                self.status_label.config(text="Configuration saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Error saving configuration: {str(e)}")
    
    def load_parameters(self):
        """Load parameters and control points from file"""
        filename = filedialog.askopenfilename(
            title="Load Configuration",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                with open(filename, 'r') as f:
                    config = json.load(f)
                
                # Update parameters
                if 'parameters' in config:
                    self.params.update(config['parameters'])
                
                # Update control points
                if config.get('slot_control_points') is not None:
                    self.slot_control_points = np.array(config['slot_control_points'])
                if config.get('flap_control_points') is not None:
                    self.flap_control_points = np.array(config['flap_control_points'])
                if config.get('hinge_point') is not None:
                    self.hinge_point = np.array(config['hinge_point'])
                if 'magnetic_snap' in config:
                    self.magnetic_snap = config['magnetic_snap']
                    if self.magnetic_snap:
                        self.magnet_style.configure("Active.TButton", background="lightgreen")
                        self.magnet_button.config(style="Active.TButton")
                    else:
                        self.magnet_button.config(style="TButton")
                
                # Update curves first
                self.slot_curve = self.bezier_curve(self.slot_control_points)
                self.flap_curve = self.bezier_curve(self.flap_control_points)
                self.assemble_components()
                
                # Update parameters from the loaded points
                self.update_parameters_from_sketch()
                
                # Plot results
                self.plot_results()
                
                self.status_label.config(text="Configuration loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Error loading configuration: {str(e)}")
    

def main():
    root = tk.Tk()
    app = AirfoilFlapDesigner(root)
    root.mainloop()

if __name__ == "__main__":
    main()