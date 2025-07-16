import numpy as np
import math
from scipy.interpolate import interp1d
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import json
import csv
import pandas as pd
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any

@dataclass
class ComponentConfig:
    """Configuration for flap and slat components"""
    flap_type: str = "Single Slotted"
    slat_type: str = "Kruger"
    has_flap: bool = True
    has_slat: bool = True

class ParameterManager:
    """Manages parameter definitions, descriptions, and validation"""
    
    def __init__(self):
        self.default_params = {
            'x1': 0.6, 'x2': 0.675, 'Delta_X1': 0.05, 'x3': 0.685, 'FractionY': 0.15,
            'DeltaX4': 0.02, 'Delta_angolo': 10, 'x6': 0.76, 'Delta_Lip': 0.045,
            'DeltaTE_main': 0.004, 'Twist Main': 0.0, 'x7': 0.84, 'DeltaX9': 0.04,
            'DeltaY9': 0.95, 'DeltaX10': 0.006, 'DeltaY10': -0.25, 'DeltaX11': 0.05,
            'x12': 0.73, 'xh_flap': 0.75, 'yh_flap': -0.1, 'Flap Deflection': 0.0,
            'x13': 0.18, 'y13': 0.08, 'x14': 0.14, 'y14': 0.06, 'x15': 0.13, 'y15': 0.03,
            'x16': 0.13, 'y16': 0.0, 'x17': 0.12, 'y17': -0.04, 'x18': 0.18, 'y18': -0.04,
            'x19': 0.14, 'y19': -0.04, 'x20': 0.11, 'y20': -0.03, 'x21': 0.1, 'y21': -0.01,
            'x22': 0.1, 'y22': 0.02, 'x23': 0.12, 'y23': 0.05, 'x24': 0.15, 'y24': 0.07,
            'xh_slat': 0.08, 'yh_slat': -0.08, 'Slat Deflection': 0.0
        }
        
        self.descriptions = {
            'x1': 'First control point X coordinate for main airfoil upper surface',
            'x2': 'Second control point X coordinate for main airfoil upper surface',
            'Delta_X1': 'Delta X adjustment for first control point',
            'x3': 'Third control point X coordinate for main airfoil',
            'FractionY': 'Fraction of Y coordinate for airfoil shaping',
            'DeltaX4': 'Delta X adjustment for fourth control point',
            'Delta_angolo': 'Delta angle adjustment for airfoil curvature',
            'x6': 'Sixth control point X coordinate',
            'Delta_Lip': 'Delta adjustment for leading edge lip',
            'DeltaTE_main': 'Delta trailing edge adjustment for main airfoil',
            'Twist Main': 'Twist angle for main airfoil section',
            'x7': 'Flap leading edge X coordinate',
            'DeltaX9': 'Delta X adjustment for flap control point 9',
            'DeltaY9': 'Delta Y adjustment for flap control point 9',
            'DeltaX10': 'Delta X adjustment for flap control point 10',
            'DeltaY10': 'Delta Y adjustment for flap control point 10',
            'DeltaX11': 'Delta X adjustment for flap control point 11',
            'x12': 'Flap trailing edge X coordinate',
            'Gap': 'Gap distance between main airfoil and flap',
            'Overlap': 'Overlap distance between main airfoil and flap',
            'xh_flap': 'Flap hinge point X coordinate',
            'yh_flap': 'Flap hinge point Y coordinate',
            'Flap Deflection': 'Flap deflection angle in degrees',
            'xh_slat': 'Slat hinge point X coordinate',
            'yh_slat': 'Slat hinge point Y coordinate',
            'Slat Deflection': 'Slat deflection angle in degrees'
        }
        
        # Add slat parameter descriptions
        for i in range(13, 25):
            self.descriptions[f'x{i}'] = f'Slat control point {i} X coordinate'
            self.descriptions[f'y{i}'] = f'Slat control point {i} Y coordinate'
        
        self.categories = {
            'Main Slot Flap - Control Points': ['x1', 'x2', 'Delta_X1', 'x3', 'FractionY', 'DeltaX4', 'Delta_angolo', 'x6', 'Delta_Lip', 'DeltaTE_main'],
            'Main Slot Slat - Control Points': [f'x{i}' for i in range(13, 19)] + [f'y{i}' for i in range(13, 19)],
            'Twist Main': ['Twist Main'],
            'Flap - Control Points': ['x7', 'DeltaX9', 'DeltaY9', 'DeltaX10', 'DeltaY10', 'DeltaX11', 'x12'],
            'Flap - Hinge & Deflection': ['Gap', 'Overlap', 'xh_flap', 'yh_flap', 'Flap Deflection'],
            'Slat - Control Points': [f'x{i}' for i in range(19, 25)] + [f'y{i}' for i in range(19, 25)],
            'Slat - Hinge & Deflection': ['xh_slat', 'yh_slat', 'Slat Deflection']
        }

class GeometryEngine:
    """Handles geometric calculations and transformations"""
    
    @staticmethod
    def bezier_curve(control_points: np.ndarray, n_points: int = 60) -> np.ndarray:
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
    
    @staticmethod
    def apply_rotation(points: np.ndarray, angle: float, center: np.ndarray) -> np.ndarray:
        """Apply rotation around a center point"""
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        translated = points - center
        rotated = np.column_stack([
            translated[:, 0] * cos_a - translated[:, 1] * sin_a,
            translated[:, 0] * sin_a + translated[:, 1] * cos_a
        ])
        return rotated + center
    
    @staticmethod
    def resample_curve(curve: np.ndarray, n_points: int) -> np.ndarray:
        """Resample curve to specified number of points"""
        if curve is None or len(curve) == 0:
            return np.array([])
        
        # Find trailing edge and reorder
        te_idx = np.argmax(curve[:, 0])
        reordered = np.concatenate([curve[te_idx:], curve[:te_idx+1]])
        
        # Create parameter for interpolation
        distances = np.concatenate([[0], np.cumsum(np.sqrt(np.diff(reordered[:, 0])**2 + np.diff(reordered[:, 1])**2))])
        t_original = distances / distances[-1]
        t_new = np.linspace(0, 1, n_points)
        
        # Interpolate
        fx = interp1d(t_original, reordered[:, 0], kind='linear')
        fy = interp1d(t_original, reordered[:, 1], kind='linear')
        
        return np.column_stack([fx(t_new), fy(t_new)])

class AirfoilDesigner:
    """Main application class"""
    
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.withdraw()
        
        # Initialize managers
        self.param_manager = ParameterManager()
        self.geometry_engine = GeometryEngine()
        
        # Initialize data structures
        self.config = ComponentConfig()
        self.params = self.param_manager.default_params.copy()
        self.airfoil_data = None
        self.f_upper = None
        self.f_lower = None
        
        # Control points and curves
        self.control_points = {}
        self.curves = {}
        self.components = {}
        
        # UI state
        self.param_vars = {}
        self.tree_item_to_param = {}
        self.dragging_point = None
        self.drag_offset = None
        self.magnetic_snap = False
        self.snap_tolerance = 0.0025
        
        # Gap/overlap tracking
        self.gap_value = 0.0
        self.overlap_value = 0.0
        self.gap_overlap_vars = {}
        
        self.setup_startup_window()
    
    def setup_startup_window(self):
        """Create and configure the startup window"""
        self.starting_window = tk.Toplevel(self.root)
        self.starting_window.title("Interactive Airfoil-Flap-Slat Designer")
        self.starting_window.geometry("350x250")
        self.starting_window.resizable(False, False)
        self.starting_window.protocol("WM_DELETE_WINDOW", self.quit_program)
        
        frame = ttk.Frame(self.starting_window, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Component selection
        self.setup_component_selection(frame)
        
        # Start button
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', pady=(30, 0))
        ttk.Button(btn_frame, text="Start Design", command=self.start_design).pack(side='right')
    
    def setup_component_selection(self, parent):
        """Setup component selection interface"""
        self.flap_type = tk.StringVar(value="Fowler")
        self.slat_type = tk.StringVar(value="Kruger")
        
        flap_types = ["None", "Fowler"]
        slat_types = ["None", "Kruger"]
        
        # Flap selection
        ttk.Label(parent, text="Flap Type:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(20, 5))
        flap_combo = ttk.Combobox(parent, width=30, values=flap_types, textvariable=self.flap_type, state='readonly')
        flap_combo.pack(pady=(0, 30), anchor='w')
        
        # Slat selection
        ttk.Label(parent, text="Slat Type:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(0, 5))
        slat_combo = ttk.Combobox(parent, width=30, values=slat_types, textvariable=self.slat_type, state='readonly')
        slat_combo.pack(pady=(0, 10), anchor='w')
    
    def start_design(self):
        """Start the main design interface"""
        # Update configuration
        self.config.has_flap = self.flap_type.get() != "None"
        self.config.has_slat = self.slat_type.get() != "None"
        self.config.flap_type = self.flap_type.get()
        self.config.slat_type = self.slat_type.get()
        
        # Setup main window
        self.starting_window.destroy()
        self.root.deiconify()
        self.root.title("Interactive Airfoil-Flap-Slat Designer")
        self.root.state('zoomed')
        
        self.setup_main_interface()
    
    def setup_main_interface(self):
        """Setup the main design interface"""
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        paned_window = ttk.PanedWindow(main_frame, orient=tk.HORIZONTAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        
        # Left panel for parameters
        left_frame = ttk.Frame(paned_window)
        paned_window.add(left_frame)
        
        # Right panel for plot
        right_frame = ttk.Frame(paned_window)
        paned_window.add(right_frame)
        
        self.setup_parameter_panel(left_frame)
        self.setup_plot_panel(right_frame)
    
    def setup_parameter_panel(self, parent):
        """Setup parameter control panel"""
        ttk.Label(parent, text="Parameters", font=('Arial', 12, 'bold')).pack(pady=(0, 10))
        
        main_container = ttk.Frame(parent)
        main_container.pack(fill=tk.BOTH, expand=True)
        
        paned_window = ttk.PanedWindow(main_container, orient=tk.VERTICAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        
        # Parameter tree
        tree_frame = ttk.Frame(paned_window)
        paned_window.add(tree_frame, weight=3)
        
        self.setup_parameter_tree(tree_frame)
        
        # Parameter details
        details_frame = ttk.Frame(paned_window)
        paned_window.add(details_frame, weight=1)
        
        self.setup_parameter_details(details_frame)
    
    def setup_parameter_tree(self, parent):
        """Setup parameter tree view"""
        tree_scroll = ttk.Scrollbar(parent)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.param_tree = ttk.Treeview(parent, yscrollcommand=tree_scroll.set, selectmode='extended')
        self.param_tree.pack(fill=tk.BOTH, expand=True)
        tree_scroll.config(command=self.param_tree.yview)
        
        # Configure columns
        self.param_tree['columns'] = ('Value',)
        self.param_tree.column('#0', width=200, minwidth=150)
        self.param_tree.column('Value', width=100, minwidth=80)
        self.param_tree.heading('#0', text='Parameter', anchor=tk.W)
        self.param_tree.heading('Value', text='Value', anchor=tk.CENTER)
        
        # Bind events
        self.param_tree.bind('<Double-1>', self.on_tree_double_click)
        self.param_tree.bind('<<TreeviewSelect>>', self.on_tree_select)
        
        self.populate_parameter_tree()
    
    def populate_parameter_tree(self):
        """Populate the parameter tree with current configuration"""
        # Filter categories based on configuration
        active_categories = self.param_manager.categories.copy()
        
        if not self.config.has_flap:
            active_categories.pop('Main Slot Flap - Control Points', None)
            active_categories.pop('Flap - Control Points', None)
            active_categories.pop('Flap - Hinge & Deflection', None)
        
        if not self.config.has_slat:
            active_categories.pop('Main Slot Slat - Control Points', None)
            active_categories.pop('Slat - Control Points', None)
            active_categories.pop('Slat - Hinge & Deflection', None)
        
        # Populate tree
        for category, params in active_categories.items():
            category_node = self.param_tree.insert('', 'end', text=category, values=('',))
            
            for param in params:
                if param in self.params:
                    value = self.params[param]
                    self.param_vars[param] = tk.DoubleVar(value=value)
                    param_node = self.param_tree.insert(category_node, 'end', text=param, 
                                                      values=(f'{value:.6f}',))
                    self.tree_item_to_param[param_node] = param
        
        # Collapse all categories
        for item in self.param_tree.get_children():
            self.param_tree.item(item, open=False)
    
    def setup_parameter_details(self, parent):
        """Setup parameter details and editor"""
        if self.config.has_flap:
            self.setup_gap_overlap_display(parent)
        
        # Update button
        ttk.Button(parent, text="Update", command=self.update_parameters).pack(fill=tk.X, pady=(20, 0))
        
        # Parameter description
        desc_frame = ttk.LabelFrame(parent, text="Parameter Details", padding=5)
        desc_frame.pack(fill=tk.X, side=tk.BOTTOM, pady=(0, 5))
        
        self.param_name_label = ttk.Label(desc_frame, text="Select a parameter", font=('Arial', 10, 'bold'))
        self.param_name_label.pack(anchor=tk.W)
        
        self.param_desc_text = tk.Text(desc_frame, height=10, width=40, wrap=tk.WORD, state='disabled')
        self.param_desc_text.pack(fill=tk.X, pady=(5, 0))
    
    def setup_gap_overlap_display(self, parent):
        """Setup gap and overlap display"""
        gap_overlap_frame = ttk.LabelFrame(parent, text="Gap & Overlap (Read-only)", padding=5)
        gap_overlap_frame.pack(fill=tk.X, pady=(10, 0))
        
        # Gap display
        gap_frame = ttk.Frame(gap_overlap_frame)
        gap_frame.pack(fill=tk.X, pady=2)
        ttk.Label(gap_frame, text="Gap:", width=8).pack(side=tk.LEFT)
        gap_var = tk.StringVar(value="0.000000")
        self.gap_overlap_vars['gap'] = gap_var
        ttk.Entry(gap_frame, textvariable=gap_var, width=12, state='readonly').pack(side=tk.LEFT, padx=(5, 0))
        
        # Overlap display
        overlap_frame = ttk.Frame(gap_overlap_frame)
        overlap_frame.pack(fill=tk.X, pady=2)
        ttk.Label(overlap_frame, text="Overlap:", width=8).pack(side=tk.LEFT)
        overlap_var = tk.StringVar(value="0.000000")
        self.gap_overlap_vars['overlap'] = overlap_var
        ttk.Entry(overlap_frame, textvariable=overlap_var, width=12, state='readonly').pack(side=tk.LEFT, padx=(5, 0))
    
    def setup_plot_panel(self, parent):
        """Setup the plot panel"""
        # Control buttons
        control_frame = ttk.Frame(parent)
        control_frame.pack(fill=tk.X, pady=(0, 5))
        
        buttons = [
            ("Restart App", self.restart_app),
            ("Load Airfoil", self.load_airfoil),
            ("Reset View", self.reset_view),
            ("Save Config", self.save_config),
            ("Load Config", self.load_config),
            ("Export", self.export_data)
        ]
        
        for text, command in buttons:
            ttk.Button(control_frame, text=text, command=command).pack(side=tk.LEFT, padx=2)
        
        # Magnetic snap button
        self.magnet_button = ttk.Button(control_frame, text="🧲 Snap", command=self.toggle_magnetic_snap)
        self.magnet_button.pack(side=tk.LEFT, padx=2)
        
        # Status label
        self.status_label = ttk.Label(control_frame, text="Load an airfoil to begin")
        self.status_label.pack(side=tk.RIGHT, padx=10)
        
        # Plot area
        self.setup_plot_area(parent)
    
    def setup_plot_area(self, parent):
        """Setup the matplotlib plot area"""
        plot_frame = ttk.Frame(parent)
        plot_frame.pack(fill=tk.BOTH, expand=True)
        
        self.fig = Figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Bind mouse events
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect("scroll_event", self.on_mouse_scroll)
    
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
                self.generate_components()
                self.plot_results()
                self.update_status()
            except Exception as e:
                messagebox.showerror("Error", f"Error loading airfoil: {str(e)}")
    
    def split_airfoil(self):
        """Split airfoil into upper and lower surfaces"""
        if self.airfoil_data is None:
            return
        
        x_le_idx = np.argmin(self.airfoil_data[:, 0])
        x_upper = self.airfoil_data[:x_le_idx+1, 0]
        y_upper = self.airfoil_data[:x_le_idx+1, 1]
        x_lower = self.airfoil_data[x_le_idx:, 0]
        y_lower = self.airfoil_data[x_le_idx:, 1]
        
        self.f_upper = interp1d(x_upper, y_upper, kind='linear', bounds_error=False, fill_value='extrapolate')
        self.f_lower = interp1d(x_lower, y_lower, kind='linear', bounds_error=False, fill_value='extrapolate')
    
    def generate_components(self):
        """Generate all components based on current configuration"""
        if self.airfoil_data is None:
            return
        
        try:
            if self.config.has_flap:
                self.generate_flap()
            if self.config.has_slat:
                self.generate_slat()
            self.assemble_components()
            if self.config.has_flap:
                self.calculate_gap_overlap()
        except Exception as e:
            self.status_label.config(text=f"Error generating components: {str(e)}")
    
    def generate_flap(self):
        """Generate flap components"""
        self.control_points['main_slot_flap'] = self.generate_main_slot_flap_control_points()
        self.control_points['flap'] = self.generate_flap_control_points()
        self.control_points['hinge_flap'] = np.array([self.params['xh_flap'], self.params['yh_flap']])
        
        self.curves['main_slot_flap'] = self.geometry_engine.bezier_curve(self.control_points['main_slot_flap'])
        self.curves['flap'] = self.geometry_engine.bezier_curve(self.control_points['flap'])
    
    def generate_slat(self):
        """Generate slat components"""
        self.control_points['main_slot_slat'] = self.generate_main_slot_slat_control_points()
        self.control_points['slat'] = self.generate_slat_control_points()
        self.control_points['hinge_slat'] = np.array([self.params['xh_slat'], self.params['yh_slat']])
        
        self.curves['main_slot_slat'] = self.geometry_engine.bezier_curve(self.control_points['main_slot_slat'])
        self.curves['slat'] = self.geometry_engine.bezier_curve(self.control_points['slat'])
    
    def generate_main_slot_flap_control_points(self) -> np.ndarray:
        """Generate main slot flap control points"""
        p = self.params
        
        # Point calculations (condensed from original)
        x1, y1 = p['x1'], self.f_lower(p['x1'])
        x2p = x1 - p['Delta_X1']
        y2p = self.f_lower(x2p)
        slope_lower = np.arctan((y2p - y1) / (x2p - x1))
        x2 = p['x2']
        y2 = y1 + slope_lower * (x2 - x1)
        
        x3 = p['x3']
        y3p = self.f_upper(x3)
        y3 = y2 + p['FractionY'] * (y3p - y2)
        
        x4 = x3 - p['DeltaX4']
        angolo_3 = 90 - 57.3 * np.arctan((y3 - y2) / (x3 - x2))
        angolo_tot = angolo_3 - p['Delta_angolo']
        slope_4 = np.tan(np.deg2rad(angolo_tot))
        y4 = y3 + slope_4 * (x4 - x3)
        
        x6 = p['x6']
        y6_upper = self.f_upper(x6)
        y6 = y6_upper - p['DeltaTE_main']
        
        x6p = x6 + 0.05
        y6p = self.f_upper(x6p)
        slope_p6 = -np.arctan((y6p - y6) / (x6p - x6))
        x5 = x6 - p['Delta_Lip']
        y5 = y6 + slope_p6 * (x6 - x5)
        
        return np.array([[x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5], [x6, y6]])
    
    def generate_flap_control_points(self) -> np.ndarray:
        """Generate flap control points"""
        p = self.params
        
        x7, y7 = p['x7'], self.f_upper(p['x7'])
        x8, y8 = p['x6'], self.f_upper(p['x6'])
        
        x9 = self.control_points['main_slot_flap'][2, 0] + p['DeltaX9']
        vdist_9 = abs(self.f_upper(p['x6']) - self.f_lower(p['x6']))
        y9 = self.control_points['main_slot_flap'][2, 1] + p['DeltaY9'] * vdist_9
        
        x10 = x9 + p['DeltaX10']
        y10 = y9 + p['DeltaY10'] * vdist_9
        
        x12, y12 = p['x12'], self.f_lower(p['x12'])
        x11 = x12 - p['DeltaX11']
        x12_bis = x12 + 0.00001
        y12_bis = self.f_lower(x12_bis)
        slope_y11 = np.arctan((y12_bis - y12) / (x12_bis - x12))
        y11 = y12 - slope_y11 * (x12 - x11)
        
        return np.array([[x7, y7], [x8, y8], [x9, y9], [x10, y10], [x11, y11], [x12, y12]])
    
    def generate_main_slot_slat_control_points(self) -> np.ndarray:
        """Generate main slot slat control points"""
        p = self.params
        return np.array([[p['x13'], p['y13']], [p['x14'], p['y14']], [p['x15'], p['y15']],
                        [p['x16'], p['y16']], [p['x17'], p['y17']], [p['x18'], p['y18']]])
    
    def generate_slat_control_points(self) -> np.ndarray:
        """Generate slat control points"""
        p = self.params
        return np.array([[p['x19'], p['y19']], [p['x20'], p['y20']], [p['x21'], p['y21']],
                        [p['x22'], p['y22']], [p['x23'], p['y23']], [p['x24'], p['y24']]])
    
    def assemble_components(self):
        """Assemble main and component surfaces"""
        if self.airfoil_data is None:
            return
        
        p = self.params
        cos_twist = np.cos(-np.deg2rad(p['Twist Main']))
        sin_twist = np.sin(-np.deg2rad(p['Twist Main']))
        
        # Get airfoil surface arrays
        x_upper = self.airfoil_data[:np.argmin(self.airfoil_data[:, 0])+1, 0]
        y_upper = self.airfoil_data[:np.argmin(self.airfoil_data[:, 0])+1, 1]
        x_lower = self.airfoil_data[np.argmin(self.airfoil_data[:, 0]):, 0]
        y_lower = self.airfoil_data[np.argmin(self.airfoil_data[:, 0]):, 1]
        
        # Assemble based on configuration
        if self.config.has_flap and self.config.has_slat:
            self.assemble_flap_slat_config(x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist)
        elif self.config.has_flap:
            self.assemble_flap_only_config(x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist)
        elif self.config.has_slat:
            self.assemble_slat_only_config(x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist)
        else:
            self.assemble_basic_config(x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist)
    
    def assemble_flap_slat_config(self, x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist):
        """Assemble configuration with both flap and slat"""
        p = self.params
        
        # Find indices
        ix6 = np.argmin(np.abs(p['x6'] - x_upper))
        ix1 = np.argmin(np.abs(p['x1'] - x_lower))
        ix7 = np.argmin(np.abs(p['x7'] - x_upper))
        ix12 = np.argmin(np.abs(p['x12'] - x_lower))
        ix13 = np.argmin(np.abs(p['x13'] - x_upper))
        ix18 = np.argmin(np.abs(p['x18'] - x_lower))
        ix19 = np.argmin(np.abs(p['x19'] - x_lower))
        ix24 = np.argmin(np.abs(p['x24'] - x_upper))
        
        # Main component
        main_x = np.concatenate([x_upper[ix6:ix13], self.curves['main_slot_slat'][:, 0], 
                                x_lower[ix18+1:ix1], self.curves['main_slot_flap'][:, 0]])
        main_y = np.concatenate([y_upper[ix6:ix13], self.curves['main_slot_slat'][:, 1], 
                                y_lower[ix18+1:ix1], self.curves['main_slot_flap'][:, 1]])
        
        self.components['main'] = self.apply_twist(main_x, main_y, cos_twist, sin_twist)
        
        # Flap component
        flap_x = np.concatenate([x_upper[:ix7], self.curves['flap'][:, 0], x_lower[ix12+1:]])
        flap_y = np.concatenate([y_upper[:ix7], self.curves['flap'][:, 1], y_lower[ix12+1:]])
        self.components['flap'] = self.apply_twist(flap_x, flap_y, cos_twist, sin_twist)
        
        # Slat component
        slat_x = np.concatenate([x_upper[ix24:-1], x_lower[:ix19], self.curves['slat'][:, 0]])
        slat_y = np.concatenate([y_upper[ix24:-1], y_lower[:ix19], self.curves['slat'][:, 1]])
        slat_x[0] = p['x24']
        slat_y[0] = self.f_upper(p['x24'])
        self.components['slat'] = self.apply_twist(slat_x, slat_y, cos_twist, sin_twist)
    
    def assemble_flap_only_config(self, x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist):
        """Assemble configuration with flap only"""
        p = self.params
        
        ix6 = np.argmin(np.abs(p['x6'] - x_upper))
        ix1 = np.argmin(np.abs(p['x1'] - x_lower))
        ix7 = np.argmin(np.abs(p['x7'] - x_upper))
        ix12 = np.argmin(np.abs(p['x12'] - x_lower))
        
        # Main component
        main_x = np.concatenate([x_upper[ix6:], x_lower[:ix1], self.curves['main_slot_flap'][:, 0]])
        main_y = np.concatenate([y_upper[ix6:], y_lower[:ix1], self.curves['main_slot_flap'][:, 1]])
        self.components['main'] = self.apply_twist(main_x, main_y, cos_twist, sin_twist)
        
        # Flap component
        flap_x = np.concatenate([x_upper[:ix7], self.curves['flap'][:, 0], x_lower[ix12+1:]])
        flap_y = np.concatenate([y_upper[:ix7], self.curves['flap'][:, 1], y_lower[ix12+1:]])
        self.components['flap'] = self.apply_twist(flap_x, flap_y, cos_twist, sin_twist)
        
        self.components['slat'] = None
    
    def assemble_slat_only_config(self, x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist):
        """Assemble configuration with slat only"""
        p = self.params
        
        ix13 = np.argmin(np.abs(p['x13'] - x_upper))
        ix18 = np.argmin(np.abs(p['x18'] - x_lower))
        ix19 = np.argmin(np.abs(p['x19'] - x_lower))
        ix24 = np.argmin(np.abs(p['x24'] - x_upper))
        
        # Main component
        main_x = np.concatenate([x_upper[:ix13], self.curves['main_slot_slat'][:, 0], x_lower[ix18+1:]])
        main_y = np.concatenate([y_upper[:ix13], self.curves['main_slot_slat'][:, 1], y_lower[ix18+1:]])
        self.components['main'] = self.apply_twist(main_x, main_y, cos_twist, sin_twist)
        
        # Slat component
        slat_x = np.concatenate([x_upper[ix24:-1], x_lower[:ix19], self.curves['slat'][:, 0]])
        slat_y = np.concatenate([y_upper[ix24:-1], y_lower[:ix19], self.curves['slat'][:, 1]])
        slat_x[0] = p['x24']
        slat_y[0] = self.f_upper(p['x24'])
        self.components['slat'] = self.apply_twist(slat_x, slat_y, cos_twist, sin_twist)
        
        self.components['flap'] = None
    
    def assemble_basic_config(self, x_upper, y_upper, x_lower, y_lower, cos_twist, sin_twist):
        """Assemble basic configuration (no flap or slat)"""
        main_x = np.concatenate([x_upper, x_lower])
        main_y = np.concatenate([y_upper, y_lower])
        self.components['main'] = self.apply_twist(main_x, main_y, cos_twist, sin_twist)
        self.components['flap'] = None
        self.components['slat'] = None
    
    def apply_twist(self, x, y, cos_twist, sin_twist):
        """Apply twist transformation to coordinates"""
        return np.column_stack([x * cos_twist - y * sin_twist, x * sin_twist + y * cos_twist])
    
    def apply_deflection(self, component_name: str) -> Optional[np.ndarray]:
        """Apply deflection to a component"""
        if component_name == 'flap' and self.components['flap'] is not None:
            angle = np.deg2rad(self.params['Flap Deflection'])
            center = self.control_points['hinge_flap']
            return self.geometry_engine.apply_rotation(self.components['flap'], -angle, center)
        
        elif component_name == 'slat' and self.components['slat'] is not None:
            angle = np.deg2rad(self.params['Slat Deflection'])
            center = self.control_points['hinge_slat']
            return self.geometry_engine.apply_rotation(self.components['slat'], angle, center)
        
        return None
    
    def calculate_gap_overlap(self):
        """Calculate gap and overlap between main and flap components"""
        if self.components['main'] is None or self.components['flap'] is None:
            self.gap_value = self.overlap_value = 0.0
            return
        
        try:
            deflected_flap = self.apply_deflection('flap')
            if deflected_flap is None:
                deflected_flap = self.components['flap']
            
            # Find main trailing edge and flap leading edge
            main_te_idx = np.argmax(self.components['main'][:, 0])
            main_te_point = self.components['main'][main_te_idx]
            flap_le_idx = np.argmin(deflected_flap[:, 0])
            
            # Calculate overlap
            x_main_max = np.max(self.components['main'][:, 0])
            x_flap_min = np.min(deflected_flap[:, 0])
            self.overlap_value = x_main_max - x_flap_min
            
            # Calculate gap (minimum distance)
            distances = np.sqrt(np.sum((main_te_point - deflected_flap)**2, axis=1))
            self.gap_value = np.min(distances)
            
            # Check for interference
            if np.any((deflected_flap[:, 0] < self.components['main'][:, 0][:, None]) & 
                     (deflected_flap[:, 1] > self.components['main'][:, 1][:, None])):
                self.gap_value = -1.0
            
        except Exception as e:
            print(f"Error calculating gap/overlap: {e}")
            self.gap_value = self.overlap_value = 0.0
        
        # Update GUI
        if hasattr(self, 'gap_overlap_vars'):
            if 'gap' in self.gap_overlap_vars:
                gap_text = "Interference" if self.gap_value < 0 else f"{self.gap_value:.6f}"
                self.gap_overlap_vars['gap'].set(gap_text)
            if 'overlap' in self.gap_overlap_vars:
                self.gap_overlap_vars['overlap'].set(f"{self.overlap_value:.6f}")
    
    def plot_results(self):
        """Plot all results"""
        if self.airfoil_data is None:
            return
        
        self.ax.clear()
        
        # Plot original airfoil
        self.ax.plot(self.airfoil_data[:, 0], self.airfoil_data[:, 1], 'k-', linewidth=1, alpha=0.5, label='Original')
        
        # Plot main component
        if self.components.get('main') is not None:
            self.ax.plot(self.components['main'][:, 0], self.components['main'][:, 1], 'b-', linewidth=2.5, label='Main')
        
        # Plot flap component
        if self.components.get('flap') is not None:
            self.plot_component_with_deflection('flap', 'r', 'Flap')
        
        # Plot slat component
        if self.components.get('slat') is not None:
            self.plot_component_with_deflection('slat', 'g', 'Slat')
        
        self.create_interactive_elements()
        self.setup_plot_appearance()
        self.canvas.draw()
    
    def plot_component_with_deflection(self, component_name: str, color: str, label: str):
        """Plot a component with optional deflection"""
        component = self.components[component_name]
        deflection_key = f'{component_name.title()} Deflection'
        deflection = self.params.get(deflection_key, 0.0)
        
        if abs(deflection) > 1e-6:
            # Plot both undeflected and deflected
            self.ax.plot(component[:, 0], component[:, 1], f'{color}--', linewidth=1.5, alpha=0.5, 
                        label=f'{label} (undeflected)')
            
            deflected = self.apply_deflection(component_name)
            if deflected is not None:
                self.ax.plot(deflected[:, 0], deflected[:, 1], f'{color}-', linewidth=2.5, 
                            label=f'{label} (deflected {deflection:.1f}°)')
        else:
            # Plot undeflected only
            self.ax.plot(component[:, 0], component[:, 1], f'{color}-', linewidth=1.5, label=f'{label} (undeflected)')
    
    def create_interactive_elements(self):
        """Create interactive control points and curves"""
        self.interactive_artists = []
        
        # Plot control points and curves for each active component
        if self.config.has_flap:
            self.plot_control_elements('main_slot_flap', 'blue', 'o')
            self.plot_control_elements('flap', 'red', 's')
            self.plot_hinge_point('hinge_flap', 'green', 'D')
        
        if self.config.has_slat:
            self.plot_control_elements('main_slot_slat', 'blue', 'o')
            self.plot_control_elements('slat', 'green', 's')
            self.plot_hinge_point('hinge_slat', 'orange', 'D')
    
    def plot_control_elements(self, element_name: str, color: str, marker: str):
        """Plot control points and connecting lines for an element"""
        if element_name not in self.control_points:
            return
        
        points = self.control_points[element_name]
        
        # Define color codes for matplotlib format strings
        color_codes = {
            'blue': 'b',
            'red': 'r',
            'green': 'g',
            'orange': 'orange',
            'purple': 'purple',
            'cyan': 'c'
        }
        
        # Define edge colors mapping
        edge_colors = {
            'blue': 'navy',
            'red': 'darkred',
            'green': 'darkgreen',
            'orange': 'darkorange',
            'purple': 'darkviolet',
            'cyan': 'darkcyan'
        }
        
        color_code = color_codes.get(color, 'k')  # Default to black
        edge_color = edge_colors.get(color, 'black')
        
        # Plot connecting lines
        for i in range(len(points) - 1):
            self.ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 
                        f'{color_code}--', linewidth=0.8, alpha=0.6, zorder=5)
        
        # Plot control points
        for point in points:
            artist = self.ax.scatter(point[0], point[1], c=color, s=80, marker=marker, 
                                   edgecolor=edge_color, linewidth=2, zorder=15, alpha=0.9)
            self.interactive_artists.append((artist, element_name))
    
    def plot_hinge_point(self, hinge_name: str, color: str, marker: str):
        """Plot a hinge point"""
        if hinge_name not in self.control_points:
            return
        
        # Define edge colors mapping
        edge_colors = {
            'blue': 'navy',
            'red': 'darkred',
            'green': 'darkgreen',
            'orange': 'darkorange',
            'purple': 'darkviolet',
            'cyan': 'darkcyan'
        }
        
        edge_color = edge_colors.get(color, 'black')
        
        point = self.control_points[hinge_name]
        artist = self.ax.scatter(point[0], point[1], c=color, s=100, marker=marker, 
                               edgecolor=edge_color, linewidth=2, zorder=15, alpha=0.9)
        self.interactive_artists.append((artist, hinge_name))
    
    def setup_plot_appearance(self):
        """Setup plot appearance and labels"""
        self.ax.set_xlabel('x/c')
        self.ax.set_ylabel('y/c')
        self.ax.grid(True, alpha=0.3)
        self.ax.axis('equal')
        self.ax.legend()
        self.ax.set_title('Interactive Airfoil-Flap-Slat Designer')
    
    def update_status(self):
        """Update status message based on current configuration"""
        if self.config.has_flap and self.config.has_slat:
            self.status_label.config(text="Flap and Slat generated! Drag control points to modify")
        elif self.config.has_flap:
            self.status_label.config(text="Flap generated! Drag control points to modify")
        elif self.config.has_slat:
            self.status_label.config(text="Slat generated! Drag control points to modify")
        else:
            self.status_label.config(text="Basic airfoil loaded successfully!")
    
    # Mouse event handlers (simplified)
    def on_mouse_press(self, event):
        """Handle mouse press events"""
        if event.inaxes != self.ax or event.button != 1:
            return
        
        # Find nearest control point
        tolerance = 0.008
        best_match = None
        min_distance = float('inf')
        
        for element_name, points in self.control_points.items():
            if isinstance(points, np.ndarray) and points.ndim == 2:
                for i, point in enumerate(points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance and distance < min_distance:
                        min_distance = distance
                        best_match = (element_name, i, point)
            elif isinstance(points, np.ndarray) and points.ndim == 1:
                distance = np.sqrt((event.xdata - points[0])**2 + (event.ydata - points[1])**2)
                if distance < tolerance and distance < min_distance:
                    min_distance = distance
                    best_match = (element_name, 0, points)
        
        if best_match:
            element_name, point_idx, point = best_match
            self.dragging_point = (element_name, point_idx)
            self.drag_offset = (event.xdata - point[0], event.ydata - point[1])
            self.status_label.config(text=f"Dragging {element_name} point {point_idx+1}")
    
    def on_mouse_move(self, event):
        """Handle mouse move events"""
        if self.dragging_point is None or event.inaxes != self.ax:
            return
        
        element_name, point_idx = self.dragging_point
        new_x = event.xdata - self.drag_offset[0]
        new_y = event.ydata - self.drag_offset[1]
        
        # Apply magnetic snap if enabled
        if self.magnetic_snap:
            new_x, new_y = self.find_nearest_airfoil_point(new_x, new_y)
        
        # Update control point
        if element_name in self.control_points:
            if isinstance(self.control_points[element_name], np.ndarray):
                if self.control_points[element_name].ndim == 2:
                    self.control_points[element_name][point_idx] = [new_x, new_y]
                else:
                    self.control_points[element_name] = np.array([new_x, new_y])
        
        self.update_display()
    
    def on_mouse_release(self, event):
        """Handle mouse release events"""
        if self.dragging_point is not None:
            self.status_label.config(text="Drag control points to modify shape")
            self.update_parameters_from_control_points()
        
        self.dragging_point = None
        self.drag_offset = None
    
    def on_mouse_scroll(self, event):
        """Handle mouse scroll for zooming"""
        if event.inaxes != self.ax:
            return
        
        scale_factor = 1.1 if event.step > 0 else 0.9
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        
        new_xlim = [(x - event.xdata) * scale_factor + event.xdata for x in xlim]
        new_ylim = [(y - event.ydata) * scale_factor + event.ydata for y in ylim]
        
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        self.canvas.draw_idle()
    
    def find_nearest_airfoil_point(self, x: float, y: float) -> Tuple[float, float]:
        """Find nearest point on original airfoil for magnetic snap"""
        if self.airfoil_data is None:
            return x, y
        
        # Create high-resolution interpolation
        t_original = np.linspace(0, 1, len(self.airfoil_data))
        t_interp = np.linspace(0, 1, len(self.airfoil_data) * 50)
        
        fx = interp1d(t_original, self.airfoil_data[:, 0], kind='linear')
        fy = interp1d(t_original, self.airfoil_data[:, 1], kind='linear')
        
        x_interp = fx(t_interp)
        y_interp = fy(t_interp)
        
        # Find nearest point
        distances = np.sqrt((x_interp - x)**2 + (y_interp - y)**2)
        min_idx = np.argmin(distances)
        
        if distances[min_idx] < self.snap_tolerance:
            return x_interp[min_idx], y_interp[min_idx]
        
        return x, y
    
    def update_display(self):
        """Update display after control point changes"""
        try:
            # Regenerate curves
            for element_name in ['main_slot_flap', 'flap', 'main_slot_slat', 'slat']:
                if element_name in self.control_points:
                    self.curves[element_name] = self.geometry_engine.bezier_curve(self.control_points[element_name])
            
            # Reassemble components
            self.assemble_components()
            
            # Update plot lines for components
            self.update_component_lines()
            
            # Update interactive elements
            self.update_interactive_elements()
            
            # Update control lines
            self.update_control_lines()
            
            # Recalculate gap/overlap if needed
            if self.config.has_flap:
                self.calculate_gap_overlap()
            
            self.canvas.draw_idle()
            
        except Exception as e:
            self.status_label.config(text=f"Error updating display: {e}")
    
    def update_component_lines(self):
        """Update the main component lines on the plot"""
        # Find and update main component line
        for line in self.ax.lines:
            if line.get_label() == 'Main' and self.components.get('main') is not None:
                line.set_data(self.components['main'][:, 0], self.components['main'][:, 1])
            elif 'Flap' in line.get_label() and self.components.get('flap') is not None:
                if 'undeflected' in line.get_label():
                    line.set_data(self.components['flap'][:, 0], self.components['flap'][:, 1])
                elif 'deflected' in line.get_label():
                    deflected = self.apply_deflection('flap')
                    if deflected is not None:
                        line.set_data(deflected[:, 0], deflected[:, 1])
                        # Update label with current deflection
                        line.set_label(f'Flap (deflected {self.params["Flap Deflection"]:.1f}°)')
            elif 'Slat' in line.get_label() and self.components.get('slat') is not None:
                if 'undeflected' in line.get_label():
                    line.set_data(self.components['slat'][:, 0], self.components['slat'][:, 1])
                elif 'deflected' in line.get_label():
                    deflected = self.apply_deflection('slat')
                    if deflected is not None:
                        line.set_data(deflected[:, 0], deflected[:, 1])
                        # Update label with current deflection
                        line.set_label(f'Slat (deflected {self.params["Slat Deflection"]:.1f}°)')
        
        # Update legend
        self.ax.legend()
    
    def update_control_lines(self):
        """Update control point connecting lines"""
        # Remove old control lines
        lines_to_remove = []
        for line in self.ax.lines:
            if line.get_linestyle() == '--' and line.get_alpha() == 0.6:
                lines_to_remove.append(line)
        
        for line in lines_to_remove:
            line.remove()
        
        # Redraw control lines
        self.draw_control_lines()
    
    def draw_control_lines(self):
        """Draw control point connecting lines"""
        # Color codes for lines
        color_codes = {
            'blue': 'b',
            'red': 'r',
            'green': 'g',
            'orange': 'orange',
            'purple': 'purple',
            'cyan': 'c'
        }
        
        # Draw lines for active components
        if self.config.has_flap:
            # Main slot flap lines
            if 'main_slot_flap' in self.control_points:
                points = self.control_points['main_slot_flap']
                for i in range(len(points) - 1):
                    self.ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 
                                'b--', linewidth=0.8, alpha=0.6, zorder=5)
            
            # Flap lines
            if 'flap' in self.control_points:
                points = self.control_points['flap']
                for i in range(len(points) - 1):
                    self.ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 
                                'r--', linewidth=0.8, alpha=0.6, zorder=5)
        
        if self.config.has_slat:
            # Main slot slat lines
            if 'main_slot_slat' in self.control_points:
                points = self.control_points['main_slot_slat']
                for i in range(len(points) - 1):
                    self.ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 
                                'b--', linewidth=0.8, alpha=0.6, zorder=5)
            
            # Slat lines
            if 'slat' in self.control_points:
                points = self.control_points['slat']
                for i in range(len(points) - 1):
                    self.ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], 
                                'g--', linewidth=0.8, alpha=0.6, zorder=5)
    
    def update_interactive_elements(self):
        """Update positions of interactive elements"""
        for artist, element_name in self.interactive_artists:
            if element_name in self.control_points:
                points = self.control_points[element_name]
                if isinstance(points, np.ndarray):
                    if points.ndim == 2:
                        artist.set_offsets(points)
                    else:
                        artist.set_offsets([points])
    
    def update_parameters_from_control_points(self):
        """Update parameters based on current control point positions"""
        try:
            # This is a simplified version - in practice, you'd need to reverse-engineer
            # the parameter values from the control point positions
            for element_name, points in self.control_points.items():
                if element_name == 'hinge_flap':
                    self.params['xh_flap'] = points[0]
                    self.params['yh_flap'] = points[1]
                elif element_name == 'hinge_slat':
                    self.params['xh_slat'] = points[0]
                    self.params['yh_slat'] = points[1]
            
            self.update_parameter_tree()
            
        except Exception as e:
            self.status_label.config(text=f"Error updating parameters: {e}")
    
    def update_parameter_tree(self):
        """Update parameter tree with current values"""
        for item_id, param_name in self.tree_item_to_param.items():
            if param_name in self.params:
                value = self.params[param_name]
                if param_name in self.param_vars:
                    self.param_vars[param_name].set(value)
                self.param_tree.set(item_id, 'Value', f'{value:.6f}')
    
    # Tree event handlers
    def on_tree_double_click(self, event):
        """Handle double-click on tree items for inline editing"""
        item = self.param_tree.selection()[0] if self.param_tree.selection() else None
        if not item or item not in self.tree_item_to_param:
            return
        
        bbox = self.param_tree.bbox(item, 'Value')
        if not bbox:
            return
        
        param_name = self.tree_item_to_param[item]
        current_value = self.param_vars[param_name].get()
        
        # Create inline editor
        self.edit_var = tk.DoubleVar(value=current_value)
        self.edit_entry = ttk.Entry(self.param_tree, textvariable=self.edit_var, width=10)
        self.edit_entry.place(x=bbox[0], y=bbox[1], width=bbox[2], height=bbox[3])
        
        self.edit_entry.select_range(0, tk.END)
        self.edit_entry.focus()
        
        self.editing_item = item
        self.editing_param = param_name
        
        # Bind completion events
        self.edit_entry.bind('<Return>', self.finish_edit)
        self.edit_entry.bind('<Escape>', self.cancel_edit)
        self.edit_entry.bind('<FocusOut>', self.finish_edit)
    
    def finish_edit(self, event=None):
        """Finish inline editing"""
        if not hasattr(self, 'edit_entry'):
            return
        
        try:
            new_value = self.edit_var.get()
            
            # Update parameter
            self.param_vars[self.editing_param].set(new_value)
            self.params[self.editing_param] = new_value
            
            # Update tree display
            self.param_tree.set(self.editing_item, 'Value', f'{new_value:.6f}')
            
            # Update components
            self.generate_components()
            self.plot_results()
            
        except tk.TclError:
            pass
        
        self.cancel_edit()
    
    def cancel_edit(self, event=None):
        """Cancel inline editing"""
        if hasattr(self, 'edit_entry'):
            self.edit_entry.destroy()
            delattr(self, 'edit_entry')
            delattr(self, 'edit_var')
            delattr(self, 'editing_item')
            delattr(self, 'editing_param')
    
    def on_tree_select(self, event):
        """Handle tree selection to show parameter details"""
        selection = self.param_tree.selection()
        if not selection:
            return
        
        item = selection[0]
        
        if item in self.tree_item_to_param:
            param_name = self.tree_item_to_param[item]
            self.param_name_label.config(text=f"Parameter: {param_name}")
            
            # Update description
            self.param_desc_text.config(state='normal')
            self.param_desc_text.delete(1.0, tk.END)
            self.param_desc_text.insert(1.0, self.param_manager.descriptions.get(param_name, "No description available"))
            self.param_desc_text.config(state='disabled')
        else:
            # Category selected
            category_name = self.param_tree.item(item, 'text')
            self.param_name_label.config(text=f"Category: {category_name}")
            self.param_desc_text.config(state='normal')
            self.param_desc_text.delete(1.0, tk.END)
            self.param_desc_text.insert(1.0, f"Select a parameter from the {category_name} category to edit its value.")
            self.param_desc_text.config(state='disabled')
    
    # Utility methods
    def update_parameters(self):
        """Update all parameters and regenerate components"""
        try:
            self.generate_components()
            self.plot_results()
            self.status_label.config(text="Parameters updated successfully!")
        except Exception as e:
            self.status_label.config(text=f"Error updating parameters: {e}")
    
    def toggle_magnetic_snap(self):
        """Toggle magnetic snap functionality"""
        self.magnetic_snap = not self.magnetic_snap
        if self.magnetic_snap:
            self.magnet_button.config(text="🧲 Snap ON")
            self.status_label.config(text="Magnetic snap ON - Points will snap to original airfoil")
        else:
            self.magnet_button.config(text="🧲 Snap OFF")
            self.status_label.config(text="Magnetic snap OFF")
    
    def reset_view(self):
        """Reset plot view to fit all data"""
        self.ax.relim()
        self.ax.autoscale()
        self.canvas.draw()
    
    def restart_app(self):
        """Restart the application"""
        self.root.destroy()
        main()
    
    def quit_program(self):
        """Quit the application"""
        self.root.destroy()
    
    # Configuration save/load
    def save_config(self):
        """Save current configuration to file"""
        filename = filedialog.asksaveasfilename(
            title="Save Configuration",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                config = {
                    'parameters': self.params,
                    'control_points': {k: v.tolist() if isinstance(v, np.ndarray) else v 
                                     for k, v in self.control_points.items()},
                    'configuration': {
                        'flap_type': self.config.flap_type,
                        'slat_type': self.config.slat_type,
                        'has_flap': self.config.has_flap,
                        'has_slat': self.config.has_slat
                    },
                    'magnetic_snap': self.magnetic_snap
                }
                
                with open(filename, 'w') as f:
                    json.dump(config, f, indent=2)
                
                self.status_label.config(text="Configuration saved successfully!")
                
            except Exception as e:
                messagebox.showerror("Error", f"Error saving configuration: {e}")
    
    def load_config(self):
        """Load configuration from file"""
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
                if 'control_points' in config:
                    for key, value in config['control_points'].items():
                        if isinstance(value, list):
                            self.control_points[key] = np.array(value)
                        else:
                            self.control_points[key] = value
                
                # Update configuration
                if 'configuration' in config:
                    cfg = config['configuration']
                    self.config.flap_type = cfg.get('flap_type', 'Fowler')
                    self.config.slat_type = cfg.get('slat_type', 'Kruger')
                    self.config.has_flap = cfg.get('has_flap', True)
                    self.config.has_slat = cfg.get('has_slat', True)
                
                # Update magnetic snap
                if 'magnetic_snap' in config:
                    self.magnetic_snap = config['magnetic_snap']
                    self.toggle_magnetic_snap()
                    self.toggle_magnetic_snap()  # Reset to correct state
                
                # Regenerate curves and components
                for element_name in ['main_slot_flap', 'flap', 'main_slot_slat', 'slat']:
                    if element_name in self.control_points:
                        self.curves[element_name] = self.geometry_engine.bezier_curve(self.control_points[element_name])
                
                self.update_parameters()
                self.assemble_components()
                self.update_parameter_tree()
                self.plot_results()
                
                self.status_label.config(text="Configuration loaded successfully!")
                
            except Exception as e:
                messagebox.showerror("Error", f"Error loading configuration: {e}")
    
    # Export functionality
    def export_data(self):
        """Open export dialog"""
        if self.components.get('main') is None:
            messagebox.showwarning("Warning", "No data to export. Please generate components first!")
            return
        
        ExportDialog(self.root, self)
    
    def export_to_csv(self, filename: str, export_options: Dict[str, Any]):
        """Export components to CSV file"""
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            for component_name in ['slat', 'main', 'flap']:
                if not export_options.get(f'export_{component_name}', False):
                    continue
                
                component = self.components.get(component_name)
                if component is None:
                    continue
                
                # Apply deflection if needed
                if component_name in ['flap', 'slat']:
                    deflected = self.apply_deflection(component_name)
                    if deflected is not None:
                        component = deflected
                
                # Resample
                n_points = export_options.get(f'{component_name}_points', 100)
                resampled = self.geometry_engine.resample_curve(component, n_points)
                
                # Write to CSV
                writer.writerow([])
                writer.writerow([f'{component_name.title()} coords (x,y,z)'])
                for point in resampled:
                    writer.writerow([point[0], point[1], 0.0])
    
    def export_to_excel(self, filename: str, export_options: Dict[str, Any]):
        """Export components to Excel file"""
        data = {'x': [], 'y': [], 'z': []}
        
        for component_name in ['slat', 'main', 'flap']:
            if not export_options.get(f'export_{component_name}', False):
                continue
            
            component = self.components.get(component_name)
            if component is None:
                continue
            
            # Apply deflection if needed
            if component_name in ['flap', 'slat']:
                deflected = self.apply_deflection(component_name)
                if deflected is not None:
                    component = deflected
            
            # Resample
            n_points = export_options.get(f'{component_name}_points', 100)
            resampled = self.geometry_engine.resample_curve(component, n_points)
            
            # Add to data
            data['x'].extend(['', f'{component_name.title()} coords (x,y,z)', ''])
            data['y'].extend(['', '', ''])
            data['z'].extend(['', '', ''])
            
            for point in resampled:
                data['x'].append(point[0])
                data['y'].append(point[1])
                data['z'].append(0.0)
        
        # Create DataFrame and save
        df = pd.DataFrame.from_dict(data, orient='index').transpose()
        df.to_excel(filename, index=False, header=False)

class ExportDialog:
    """Dialog for export settings"""
    
    def __init__(self, parent, designer):
        self.parent = parent
        self.designer = designer
        
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Export Settings")
        self.dialog.geometry("300x350")
        self.dialog.resizable(False, False)
        
        self.setup_dialog()
    
    def setup_dialog(self):
        """Setup export dialog interface"""
        frame = ttk.Frame(self.dialog, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # File format selection
        ttk.Label(frame, text="Export Format:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(0, 5))
        self.file_format = tk.StringVar(value="CSV")
        format_combo = ttk.Combobox(frame, values=["CSV", "XLSX"], state='readonly', textvariable=self.file_format)
        format_combo.pack(anchor='w', pady=(0, 20))
        
        # Component selection
        ttk.Label(frame, text="Export Components:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(0, 5))
        
        self.export_vars = {}
        self.point_vars = {}
        
        components = [
            ('slat', 'Slat', 75, self.designer.config.has_slat),
            ('main', 'Main', 150, True),
            ('flap', 'Flap', 75, self.designer.config.has_flap)
        ]
        
        for comp_name, comp_label, default_points, enabled in components:
            # Checkbox
            var = tk.BooleanVar(value=enabled)
            self.export_vars[comp_name] = var
            check = ttk.Checkbutton(frame, text=f"Export {comp_label}", variable=var, state='normal' if enabled else 'disabled')
            check.pack(anchor='w', pady=2)
            
            # Points entry
            points_frame = ttk.Frame(frame)
            points_frame.pack(fill='x', pady=2)
            
            label = ttk.Label(points_frame, text=f"{comp_label} points:")
            label.pack(side='left')
            
            points_var = tk.IntVar(value=default_points)
            self.point_vars[comp_name] = points_var
            entry = ttk.Entry(points_frame, textvariable=points_var, width=8, state='normal' if enabled else 'disabled')
            entry.pack(side='right')
        
        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', pady=(30, 0))
        
        ttk.Button(btn_frame, text="Export", command=self.perform_export).pack(side='right', padx=(5, 0))
        ttk.Button(btn_frame, text="Cancel", command=self.dialog.destroy).pack(side='right')
    
    def perform_export(self):
        """Perform the export operation"""
        # Collect export options
        export_options = {}
        for comp_name, var in self.export_vars.items():
            export_options[f'export_{comp_name}'] = var.get()
        
        for comp_name, var in self.point_vars.items():
            export_options[f'{comp_name}_points'] = var.get()
        
        # Check if at least one component is selected
        if not any(export_options[f'export_{comp}'] for comp in ['slat', 'main', 'flap']):
            messagebox.showwarning("Warning", "Please select at least one component to export!")
            return
        
        # Get filename
        file_format = self.file_format.get()
        if file_format == "CSV":
            filename = filedialog.asksaveasfilename(
                title="Export to CSV",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
            )
        else:
            filename = filedialog.asksaveasfilename(
                title="Export to Excel",
                defaultextension=".xlsx",
                filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
            )
        
        if filename:
            try:
                if file_format == "CSV":
                    self.designer.export_to_csv(filename, export_options)
                else:
                    self.designer.export_to_excel(filename, export_options)
                
                self.dialog.destroy()
                self.designer.status_label.config(text="Export completed successfully!")
                
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {e}")


def main():
    """Main entry point"""
    root = tk.Tk()
    app = AirfoilDesigner(root)
    root.mainloop()


if __name__ == "__main__":
    main()