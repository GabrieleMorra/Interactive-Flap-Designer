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
import pandas as pd
import openpyxl
from functools import partial
from tkinter import simpledialog

class AirfoilFlapDesigner:

    def __init__(self, root):
        """Initialize the Airfoil Flap Designer GUI"""

        # Initialize pre-main window
        self.root = root
        self.root.withdraw()  # Hide the root window initially

        self.starting_window = tk.Toplevel(self.root)
        self.starting_window.title("Interactive Airfoil-Flap-Slat Designer")
        self.starting_window.geometry("350x250")
        self.starting_window.resizable(False, False)
        
        self.starting_window.protocol("WM_DELETE_WINDOW", self.quit_program)

        # Variables
        self.create_slat = tk.BooleanVar(value=True)
        self.create_flap = tk.BooleanVar(value=True)
        self.slat_type = tk.StringVar(value="None")
        self.flap_type = tk.StringVar(value="Fowler")

        flap_types = ["None", "Fowler"]
        slat_types = ["None", "Kruger"]
        
        # Main frame
        frame = ttk.Frame(self.starting_window, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        self.flap_type = tk.StringVar()
        self.slat_type = tk.StringVar()
        
        # Starting options
        self.flap_label = ttk.Label(frame, text="Flap Type:", font=('Arial', 10, 'bold'))
        self.flap_label.pack(anchor='w', pady=(20, 5))
        flap_selection = ttk.Combobox(frame, width=30, values=flap_types, textvariable=self.flap_type, state='readonly')
        flap_selection.set("Fowler")
        flap_selection.pack(pady=(0, 30), anchor='w')
        
        self.slat_label = ttk.Label(frame, text="Slat Type:", font=('Arial', 10, 'bold'))
        self.slat_label.pack(anchor='w', pady=(0, 5))
        slat_selection = ttk.Combobox(frame, width=30, values=slat_types, textvariable=self.slat_type, state='readonly')
        slat_selection.set("Kruger")
        slat_selection.pack(pady=(0, 10), anchor='w')

        
        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', pady=(30,0))

        ttk.Button(btn_frame, text="Start Design", command=self.start_design).pack(side='right', padx=(0,0))

    def quit_program(self):
        self.root.destroy()

    def start_design(self):
        self.starting_window.destroy()  # Close the starting window
        self.root.deiconify()  # Show the main window
        self.root.title("Interactive Airfoil-Flap-Slat Designer")
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
            'DeltaTE_main': 0.004, 'Twist Main': 0.0, 'x7': 0.84, 'DeltaX9': 0.04,
            'DeltaY9': 0.95, 'DeltaX10': 0.006, 'DeltaY10': -0.25, 'DeltaX11': 0.05,
            'x12': 0.73, 'xh_flap': 0.75, 'yh_flap': -0.1, 'Flap Deflection': 0.0,
            'x13': 0.18, 'y13': 0.08, 'x14': 0.14, 'y14': 0.06, 'x15': 0.13, 'y15': 0.03,
            'x16': 0.13, 'y16': 0.0, 'x17': 0.12, 'y17': -0.04, 'x18': 0.18, 'y18': -0.04,
            'x19': 0.14, 'y19': -0.04, 'x20': 0.11, 'y20': -0.03, 'x21': 0.1, 'y21': -0.01,
            'x22': 0.1, 'y22': 0.02, 'x23': 0.12, 'y23': 0.05, 'x24': 0.15, 'y24': 0.07,
            'xh_slat': 0.08, 'yh_slat': -0.08, 'Slat Deflection': 0.0
        }
        
        # Control points and interactive elements
        self.main_slot_flap_control_points = None
        self.flap_control_points = None
        self.flap_hinge_point = None
        self.main_slot_flap_artists = []
        self.flap_artists = []
        self.flap_hinge_artist = None

        self.point_annotations = []
        self.dragging_point = None
        self.drag_offset = None

        self.main_slot_slat_control_points = None
        self.slat_control_points = None
        self.slat_hinge_point = None
        self.main_slot_slat_artists = []
        self.slat_artists = []
        self.slat_hinge_artist = None

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
        self.main_slot_flap_curve = None
        self.main_slot_slat_curve = None
        self.flap_curve = None
        
        # Plot elements
        self.main_slot_flap_curve_line = None
        self.main_slot_slat_curve_line = None
        self.flap_curve_line = None
        self.main_line = None
        self.flap_line = None
        self.slat_line = None
        
        self.setup_main_gui()
    
    def setup_main_gui(self):
        # Create main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create PanedWindow
        paned_window = ttk.PanedWindow(main_frame, orient=tk.HORIZONTAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        
        # Left panel for parameters - FORZA la larghezza iniziale
        left_frame = ttk.Frame(paned_window)
        paned_window.add(left_frame)
        
        # Right panel for plot
        right_frame = ttk.Frame(paned_window)
        paned_window.add(right_frame)
        
        # IMPORTANTE: Assicurati che il left_frame si espanda
        left_frame.pack_propagate(True)  # Permetti espansione
        
        self.setup_parameter_panel(left_frame)
        self.setup_plot_panel(right_frame)

    def setup_parameter_panel(self, parent):
        # Title
        ttk.Label(parent, text="Parameters", font=('Arial', 12, 'bold')).pack(pady=(0, 10))
        
        # Main container - single column layout
        main_container = ttk.Frame(parent)
        main_container.pack(fill=tk.BOTH, expand=True)

        paned_window = ttk.PanedWindow(main_container, orient=tk.VERTICAL)
        paned_window.pack(fill=tk.BOTH, expand=True)
        
        # Top half - TreeView
        tree_frame = ttk.Frame(paned_window)
        paned_window.add(tree_frame, weight=3)
        
        # TreeView with scrollbar
        tree_scroll = ttk.Scrollbar(tree_frame)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.param_tree = ttk.Treeview(tree_frame, yscrollcommand=tree_scroll.set, selectmode='extended')
        self.param_tree.pack(fill=tk.BOTH, expand=True)
        tree_scroll.config(command=self.param_tree.yview)
        self.param_tree.bind('<Double-1>', self.on_tree_double_click)
        
        # Configure TreeView columns
        self.param_tree['columns'] = ('Value',)
        self.param_tree.column('#0', width=200, minwidth=150)
        self.param_tree.column('Value', width=100, minwidth=80)
        
        self.param_tree.heading('#0', text='Parameter', anchor=tk.W)
        self.param_tree.heading('Value', text='Value', anchor=tk.CENTER)
        
        # Bottom half - Parameter details and editor
        details_frame = ttk.Frame(main_container)
        paned_window.add(details_frame, weight=1)

        if self.flap_type.get() != "None":
            # Gap and Overlap display section
            gap_overlap_frame = ttk.LabelFrame(details_frame, text="Gap & Overlap (Read-only)", padding=5)
            gap_overlap_frame.pack(fill=tk.X, pady=(10, 0))
            
            # Gap display
            gap_frame = ttk.Frame(gap_overlap_frame)
            gap_frame.pack(fill=tk.X, pady=2)
            ttk.Label(gap_frame, text="Gap:", width=8).pack(side=tk.LEFT)
            gap_var = tk.StringVar(value="0.000000")
            self.gap_overlap_vars = {'gap': gap_var}
            gap_entry = ttk.Entry(gap_frame, textvariable=gap_var, width=12, state='readonly')
            gap_entry.pack(side=tk.LEFT, padx=(5, 0))
            
            # Overlap display
            overlap_frame = ttk.Frame(gap_overlap_frame)
            overlap_frame.pack(fill=tk.X, pady=2)
            ttk.Label(overlap_frame, text="Overlap:", width=8).pack(side=tk.LEFT)
            overlap_var = tk.StringVar(value="0.000000")
            self.gap_overlap_vars['overlap'] = overlap_var
            overlap_entry = ttk.Entry(overlap_frame, textvariable=overlap_var, width=12, state='readonly')
            overlap_entry.pack(side=tk.LEFT, padx=(5, 0))
        
        # Parameter description section
        desc_frame = ttk.LabelFrame(details_frame, text="Parameter Details", padding=5)
        desc_frame.pack(fill=tk.X, side=tk.BOTTOM, pady=(0, 5))
        
        # Parameter name
        self.param_name_label = ttk.Label(desc_frame, text="Select a parameter", font=('Arial', 10, 'bold'))
        self.param_name_label.pack(anchor=tk.W)
        
        # Parameter description
        self.param_desc_text = tk.Text(desc_frame, height=10, width=40, wrap=tk.WORD, state='disabled')
        self.param_desc_text.pack(fill=tk.X, pady=(5, 0))
        
        ##############################################################
        ## Older version: parameter was edited in a separate dialog ##
        # # Parameter editor section
        # editor_frame = ttk.LabelFrame(details_frame, text="Edit Parameter", padding=5)
        # editor_frame.pack(fill=tk.X, pady=(5, 0))
        
        # # Value entry
        # self.value_frame = ttk.Frame(editor_frame)
        # self.value_frame.pack(fill=tk.X, pady=2)
        
        # ttk.Label(self.value_frame, text="Value:", width=8).pack(side=tk.LEFT)
        # self.param_value_var = tk.DoubleVar()
        # self.param_value_entry = ttk.Entry(self.value_frame, textvariable=self.param_value_var, width=15)
        # self.param_value_entry.pack(side=tk.LEFT, padx=(5, 0))
        # self.param_value_entry.bind('<Return>', lambda e: self.update_parameters())
        
        # # Update All button at bottom
        # ttk.Button(details_frame, text="Update Parameters", 
        #         command=self.update_parameters).pack(fill=tk.X, pady=(10, 0))
        ##############################################################
        
        # Initialize parameter data
        self.init_parameter_tree()
        
        # Bind tree selection event
        self.param_tree.bind('<<TreeviewSelect>>', self.on_tree_select)

    def on_tree_double_click(self, event):
        """Handle double-click on tree items to edit parameter value inline"""
        item = self.param_tree.selection()[0] if self.param_tree.selection() else None
        if not item:
            return
        
        # Check if it's a parameter (not a category)
        if item not in self.tree_item_to_param:
            return
        
        # Get the bounding box of the Value column
        bbox = self.param_tree.bbox(item, 'Value')
        if not bbox:
            return
        
        param_name = self.tree_item_to_param[item]
        current_value = self.param_vars[param_name].get()
        
        # Create an Entry widget positioned over the cell
        self.edit_var = tk.DoubleVar(value=current_value)
        self.edit_entry = ttk.Entry(self.param_tree, textvariable=self.edit_var, width=10)
        
        # Position the entry over the cell
        self.edit_entry.place(x=bbox[0], y=bbox[1], width=bbox[2], height=bbox[3])
        
        # Select all text and focus
        self.edit_entry.select_range(0, tk.END)
        self.edit_entry.focus()
        
        # Store current editing info
        self.editing_item = item
        self.editing_param = param_name
        
        # Bind events to finish editing
        self.edit_entry.bind('<Return>', self.finish_edit)
        self.edit_entry.bind('<Escape>', self.cancel_edit)
        self.edit_entry.bind('<FocusOut>', self.finish_edit)

    def finish_edit(self, event=None):
        """Finish editing and update the parameter"""
        if not hasattr(self, 'edit_entry'):
            return
        
        try:
            new_value = self.edit_var.get()
            
            # Update the parameter
            self.param_vars[self.editing_param].set(new_value)
            self.params[self.editing_param] = new_value
            
            # Update tree display
            self.param_tree.set(self.editing_item, 'Value', f'{new_value:.6f}')
            
            # Update the entry field if this parameter is currently selected
            if hasattr(self, 'current_param') and self.current_param == self.editing_param:
                self.param_value_var.set(new_value)
            
            # Update the plot/calculation
            self.update_parameters()
            
        except tk.TclError:
            pass  # Invalid value, ignore
        
        # Clean up
        self.cancel_edit()

    def cancel_edit(self, event=None):
        """Cancel editing and remove the entry widget"""
        if hasattr(self, 'edit_entry'):
            self.edit_entry.destroy()
            del self.edit_entry
            del self.edit_var
            del self.editing_item
            del self.editing_param
        
    def init_parameter_tree(self):
        """Initialize the parameter tree with categories and parameters"""
        
        # Parameter descriptions
        self.param_descriptions = {
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
        
        # Add slat parameters descriptions
        for i in range(13, 25):
            self.param_descriptions[f'x{i}'] = f'Slat control point {i} X coordinate'
            self.param_descriptions[f'y{i}'] = f'Slat control point {i} Y coordinate'
        
        # Parameter categories
        tree_categories = {
            'Main Slot Flap - Control Points': ['x1', 'x2', 'Delta_X1', 'x3', 'FractionY', 'DeltaX4', 'Delta_angolo', 'x6', 'Delta_Lip', 'DeltaTE_main'],
            'Main Slot Slat - Control Points': [f'x{i}' for i in range(13, 19)] + [f'y{i}' for i in range(13, 19)],
            'Twist Main': ['Twist Main'],
            'Flap - Control Points': ['x7', 'DeltaX9', 'DeltaY9', 'DeltaX10', 'DeltaY10', 'DeltaX11', 'x12'],
            'Flap - Hinge & Defletion': ['Gap', 'Overlap', 'xh_flap', 'yh_flap', 'Flap Deflection'],
            'Slat - Control Points': [f'x{i}' for i in range(19, 25)] + [f'y{i}' for i in range(19, 25)],
            'Slat - Hinge & Defletion': ['xh_slat', 'yh_slat', 'Slat Deflection']
        }

        if self.flap_type.get() == "None":
            del tree_categories['Main Slot Flap - Control Points']
            del tree_categories['Flap - Control Points']
            del tree_categories['Flap - Hinge & Defletion']
        if self.slat_type.get() == "None":
            del tree_categories['Main Slot Slat - Control Points']
            del tree_categories['Slat - Control Points']
            del tree_categories['Slat - Hinge & Defletion']
        
        # Store parameter variables and tree item mapping
        self.param_vars = {}
        self.tree_item_to_param = {}  # Map tree items to parameter names
        
        # Populate tree
        for category, params in tree_categories.items():
            # Create category node
            category_node = self.param_tree.insert('', 'end', text=category, values=('',))
            
            # Add parameters to category
            for param in params:
                if param in self.params:
                    value = self.params[param]
                    self.param_vars[param] = tk.DoubleVar(value=value)
                    param_node = self.param_tree.insert(category_node, 'end', text=param, 
                                                    values=(f'{value:.6f}',))
                    # Store parameter name mapping
                    self.tree_item_to_param[param_node] = param
        
        # Expand all categories by default - CHANGED TO COLLAPSED
        for item in self.param_tree.get_children():
            self.param_tree.item(item, open=False)  # Keep categories collapsed

    def update_tree_parameter(self):
        """Update the parameter tree with new value for a given parameter"""
        selection = self.param_tree.selection()
        if not selection:
            return
        
        param_name = selection[0]
        new_value = self.param_value_var.get()

        for item_id in self.tree_item_to_param.items():
            if item_id[0] == param_name:
                self.param_tree.set(item_id[0], 'Value', f'{new_value:.6f}')
                break

    def update_tree_parameters(self):
        """Update all tree parameters with current values"""
        for item_id, param_name in self.tree_item_to_param.items():
            if param_name in self.param_vars:
                value = self.params[param_name]
                self.param_vars[param_name].set(value)
                self.param_tree.set(item_id, 'Value', f'{value:.6f}')
            else:
                # If parameter is not in current vars, set to empty
                self.param_tree.set(item_id, 'Value', '')

        if hasattr(self, 'current_param'):
            if self.current_param is not None:
                self.param_value_var.set(self.params[self.current_param])
            

    def on_tree_select(self, event):
        """Handle tree selection to show parameter details"""
        selection = self.param_tree.selection()
        if not selection:
            return
        
        item = selection[0]
        
        # Check if it's a parameter (has mapping in tree_item_to_param)
        if item in self.tree_item_to_param:
            param_name = self.tree_item_to_param[item]
            
            # Update parameter name label
            self.param_name_label.config(text=f"Parameter: {param_name}")
            
            # Update description
            self.param_desc_text.config(state='normal')
            self.param_desc_text.delete(1.0, tk.END)
            self.param_desc_text.insert(1.0, self.param_descriptions[param_name])
            self.param_desc_text.config(state='disabled')
            
            # Update value entry
            if param_name in self.param_vars:
                self.param_value_var.set(self.param_vars[param_name].get())
                self.param_value_entry.config(state='normal')
                self.current_param = param_name
            else:
                self.param_value_entry.config(state='disabled')
                self.current_param = None
        else:
            # Category selected
            category_name = self.param_tree.item(item, 'text')
            self.param_name_label.config(text=f"Category: {category_name}")
            self.param_desc_text.config(state='normal')
            self.param_desc_text.delete(1.0, tk.END)
            self.param_desc_text.insert(1.0, f"Select a parameter from the {category_name} category to edit its value.")
            self.param_desc_text.config(state='disabled')
            self.param_value_entry.config(state='disabled')
            self.current_param = None

    def setup_plot_panel(self, parent):
        # Top control panel
        control_frame = ttk.Frame(parent)
        control_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Button(control_frame, text="Restart App", command=self.restart_code).pack(side=tk.LEFT, padx=2)
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
        
        self.fig = Figure(figsize=(8, 8))
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
        t_interp = np.linspace(0, 1, len(x_curve) * 50)  # 50x più punti
        
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
        export_window.geometry("300x330")
        export_window.resizable(False, False)
        
        # Variables
        self.export_main = tk.BooleanVar(value=True)
        self.export_flap = tk.BooleanVar(value=True)
        self.export_slat = tk.BooleanVar(value=True)
        self.slat_points = tk.IntVar(value=75)
        self.main_points = tk.IntVar(value=150)
        self.flap_points = tk.IntVar(value=75)
        
        # Main frame
        frame = ttk.Frame(export_window, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        
        # Export options
        ttk.Label(frame, text="Export Options:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(0,5))

        self.file_type = tk.StringVar()
        file_type = ttk.Combobox(frame, values=["CSV", "XLSX"], state='readonly', textvariable=self.file_type)
        file_type.set("CSV")
        file_type.pack(anchor='w', pady=(0, 20))

        slat_check = ttk.Checkbutton(frame, text="Export Slat Component", variable=self.create_slat,
                                    command=self.update_export_options)
        slat_check.pack(anchor='w', pady=(0, 3))
        if self.slat_type.get() == "None":
            slat_check.config(state='disabled')
            self.export_slat.set(False)
        
        main_check = ttk.Checkbutton(frame, text="Export Main Component", variable=self.export_main,
                                    command=self.update_export_options)
        main_check.pack(anchor='w', pady=(0, 3))
        
        flap_check = ttk.Checkbutton(frame, text="Export Flap Component", variable=self.export_flap,
                                    command=self.update_export_options)
        flap_check.pack(anchor='w', pady=(0, 3))
        if self.flap_type.get() == "None":
            flap_check.config(state='disabled')
            self.export_flap.set(False)

        
        # Points selection
        ttk.Label(frame, text="Number of Points:", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(10,5))

        # Slat points
        slat_frame = ttk.Frame(frame)
        slat_frame.pack(fill='x', pady=2)
        self.slat_label = ttk.Label(slat_frame, text="Slat points:")
        self.slat_label.pack(side='left')
        self.slat_entry = ttk.Entry(slat_frame, textvariable=self.slat_points, width=8)
        self.slat_entry.pack(side='right')
        if self.slat_type.get() == "None":
            self.slat_entry.config(state='disabled')
            self.slat_label.config(foreground='gray')
        
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
        if self.flap_type.get() == "None":
            self.flap_entry.config(state='disabled')
            self.flap_label.config(foreground='gray')
        
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
            min_distance = float('inf')
            
            # Calcola la distanza minima tra il trailing edge e tutti i punti del flap
            for flap_point in deflected_flap:
                distance = np.sqrt((main_te_point[0] - flap_point[0])**2 + 
                                (main_te_point[1] - flap_point[1])**2)
                if distance < min_distance:
                    min_distance = distance
                # Se il punto è a x minori e y maggiori di un qualsiasi punto del main, considera gap -1
                cond = np.logical_and(
                    flap_point[0] < self.main_component[:, 0],
                    flap_point[1] > self.main_component[:, 1]
                )
                if np.any(cond):
                    min_distance = -1.0
                    break

            self.gap_value = min_distance
            
                
        except Exception as e:
            print(f"Error calculating gap/overlap: {str(e)}")
            self.gap_value = 0.0
            self.overlap_value = 0.0
        
        # Update GUI variables
        if hasattr(self, 'gap_overlap_vars'):
            if 'gap' in self.gap_overlap_vars:
                if self.gap_value < 0:
                    self.gap_overlap_vars['gap'].set("Interference")
                else:
                    self.gap_overlap_vars['gap'].set(f"{self.gap_value:.6f}")

            if 'overlap' in self.gap_overlap_vars:
                self.gap_overlap_vars['overlap'].set(f"{self.overlap_value:.6f}")

    def update_export_options(self):
        """Update export options based on checkboxes"""
        main_enabled = self.export_main.get()
        flap_enabled = self.export_flap.get()
        slat_enabled = self.export_slat.get()
        
        # Enable/disable entries
        self.main_entry.config(state='normal' if main_enabled else 'disabled')
        self.flap_entry.config(state='normal' if flap_enabled else 'disabled')
        self.slat_entry.config(state='normal' if slat_enabled else 'disabled')
        self.main_label.config(foreground='black' if main_enabled else 'gray')
        self.flap_label.config(foreground='black' if flap_enabled else 'gray')
        self.slat_label.config(foreground='black' if slat_enabled else 'gray')

    def perform_export(self, export_window):
        """Perform the actual export"""
        if not self.export_main.get() and not self.export_flap.get():
            messagebox.showwarning("Warning", "Please select at least one component to export!")
            return
    
        if self.file_type.get() == "CSV":
            filename = filedialog.asksaveasfilename(
                title="Export to CSV",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
            )
        elif self.file_type.get() == "XLSX":
            filename = filedialog.asksaveasfilename(
                title="Export to XLSX",
                defaultextension=".xlsx",
                filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
            )
        
        if filename:
            try:
                if self.file_type.get() == "XLSX":
                    self.export_to_excel(filename)

                elif self.file_type.get() == "CSV":
                    self.export_to_csv(filename)

                export_window.destroy()
                self.status_label.config(text="Export completed successfully!")

            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {str(e)}")

    def export_to_excel(self, filename):
        """Export components to Excel file"""
        data = {
            'x': [],
            'y': [],
            'z': []
        }

        if self.export_slat.get():
            slat_deflected = self.apply_slat_deflection()
            slat_data = self.resample_component(slat_deflected, self.slat_points.get())
            data['x'].append("")
            data['y'].append("")
            data['z'].append("")
            data['x'].append("Slat coords (x,y,z)")
            data['y'].append("")
            data['z'].append("")
            for i in range(slat_data.shape[0]):
                data['x'].append(slat_data[i, 0].tolist())
                data['y'].append(slat_data[i, 1].tolist())
                data['z'].append(0.0) 

        if self.export_main.get():
            main_data = self.resample_component(self.main_component, self.main_points.get())
            data['x'].append("")
            data['y'].append("")
            data['z'].append("")
            data['x'].append("Main coords (x,y,z)")
            data['y'].append("")
            data['z'].append("")
            for i in range(main_data.shape[0]):
                data['x'].append(main_data[i, 0].tolist())
                data['y'].append(main_data[i, 1].tolist())
                data['z'].append(0.0) 

        if self.export_flap.get():
            flap_deflected = self.apply_flap_deflection()
            flap_data = self.resample_component(flap_deflected, self.flap_points.get())
            data['x'].append("")
            data['y'].append("")
            data['z'].append("")
            data['x'].append("Flap coords (x,y,z)")
            data['y'].append("")
            data['z'].append("")
            for i in range(flap_data.shape[0]):
                data['x'].append(flap_data[i, 0].tolist())
                data['y'].append(flap_data[i, 1].tolist())
                data['z'].append(0.0)
        
        
        df = pd.DataFrame.from_dict(data, orient='index').transpose()
        df.to_excel(filename, index=False, header=False)

    def export_to_csv(self, filename):
        """Export components to CSV file"""
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')

            if self.export_slat.get():
                writer.writerow("")
                writer.writerow(['Slat coords (x,y)'])
                slat_deflected = self.apply_slat_deflection()
                slat_deflected = self.resample_component(slat_deflected, self.slat_points.get())
                for point in slat_deflected:
                    writer.writerow([point[0], point[1], 0.0])
            
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

    def apply_slat_deflection(self):
        """Apply rigid rotation to slat around hinge point"""

        if self.slat_component is None or abs(self.params['Slat Deflection']) < 1e-6:
            return self.slat_component.copy() if self.slat_component is not None else None
        
        xh, yh = self.params['xh_slat'], self.params['yh_slat']
        delta = np.deg2rad(self.params['Slat Deflection'])  # Positive = nose down

        # Translate slat so hinge is at origin
        x_trans = self.slat_component[:, 0] - xh
        y_trans = self.slat_component[:, 1] - yh

        # Apply rotation matrix (counterclockwise for positive deflection)
        cos_delta = np.cos(delta) # Positive for counterclockwise rotation
        sin_delta = np.sin(delta) 

        x_rot = x_trans * cos_delta - y_trans * sin_delta
        y_rot = x_trans * sin_delta + y_trans * cos_delta

        deflected_slat = np.column_stack([x_rot + xh, y_rot + yh])

        return deflected_slat

    def apply_flap_deflection(self):
        """Apply rigid rotation to flap around hinge point"""
        if self.flap_component is None or abs(self.params['Flap Deflection']) < 1e-6:
            return self.flap_component.copy() if self.flap_component is not None else None
        
        # Get hinge point and deflection angle
        xh, yh = self.params['xh_flap'], self.params['yh_flap']
        delta = np.deg2rad(self.params['Flap Deflection'])  # Positive = nose down
        
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

    def update_parameters(self):
        """Update all parameters and regenerate"""
        try:
            # Salva i limiti correnti
            current_xlim = self.ax.get_xlim()
            current_ylim = self.ax.get_ylim()
            
            self.params[self.current_param] = float(self.param_value_entry.get())

            self.update_tree_parameters()
            
            if self.flap_type.get() != "None" and self.slat_type.get() != "None":
                    self.status_label.config(text="Flap and Slat generated! Drag control points to modify")
                    self.generate_flap()
                    self.generate_slat()
                    self.assemble_components()
                    self.plot_results()
                    self.calculate_gap_overlap()
                    

            elif self.flap_type.get() != "None" and self.slat_type.get() == "None":
                self.status_label.config(text="Flap generated! Drag control points to modify")
                self.generate_flap()
                self.assemble_components()
                self.plot_results()
                self.calculate_gap_overlap()
                
                
            elif self.flap_type.get() == "None" and self.slat_type.get() != "None":
                self.status_label.config(text="Slat generated! Drag control points to modify")
                self.generate_slat()
                self.assemble_components()
                self.plot_results()
                
            
            elif self.flap_type.get() == "None" and self.slat_type.get() == "None":
                self.status_label.config(text="Basic airfoil loaded successfully!")
                self.assemble_components()
                self.plot_results()
                
            
            # Ripristina i limiti
            self.ax.set_xlim(current_xlim)
            self.ax.set_ylim(current_ylim)
            self.canvas.draw()
        except Exception as e:
            self.status_label.config(text=f"Error updating parameters: {str(e)}")


    def update_parameters_from_sketch(self):
        """Update all parameters based on control points moved in the sketch"""
        try:
            if self.main_slot_flap_control_points is None or self.flap_control_points is None:
                return
            
            p = self.params

            # Extract control points
            x1, y1 = self.main_slot_flap_control_points[0]
            x2, y2 = self.main_slot_flap_control_points[1] 
            x3, y3 = self.main_slot_flap_control_points[2]
            x4, y4 = self.main_slot_flap_control_points[3]
            x5, y5 = self.main_slot_flap_control_points[4]
            x6, y6 = self.main_slot_flap_control_points[5]
            
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
            p['xh_flap'] = self.flap_hinge_point[0]
            p['yh_flap'] = self.flap_hinge_point[1]
            
            # Delta_X1: calcolo dalla geometria attuale
            y1_orig = self.f_lower(x1)
            # Calcola la pendenza della curva originale in x1
            dx = 0.001
            y1_left = self.f_lower(x1 - dx)
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
            y3_upper = self.f_upper(x3)
            if abs(y3_upper - y2) > 1e-6:
                p['FractionY'] = max(0, min(1, (y3 - y2) / (y3_upper - y2)))
            
            # DeltaX4
            p['DeltaX4'] = abs(x3 - x4)
            
            # Delta_angolo: CORREZIONE QUI
            if x3 != x2 and x4 != x3:
                # Calcola angolo_3 come nel codice originale (generate_main_slot_flap_control_points)
                angolo_3 = 90 - 57.3 * np.arctan((y3 - y2) / (x3 - x2))
                
                # Calcola angolo_tot dal vettore 3->4
                slope_34 = (y4 - y3) / (x4 - x3)
                angolo_tot = np.degrees(np.arctan(slope_34))
                
                # Delta_angolo = angolo_3 - angolo_tot (come nella definizione originale)
                p['Delta_angolo'] = angolo_3 - angolo_tot
            
            # Delta_Lip
            p['Delta_Lip'] = abs(x6 - x5)
            
            # DeltaTE_main
            y6_upper = self.f_upper(x6)
            p['DeltaTE_main'] = abs(y6_upper - y6)
            
            # Flap parameters
            p['DeltaX9'] = x9 - x3
            
            # DeltaY9: normalizzato rispetto alla distanza verticale
            vdist = abs(self.f_upper(x6) - self.f_lower(x6))
            if vdist > 1e-6:
                p['DeltaY9'] = (y9 - y3) / vdist
            
            # DeltaX10
            p['DeltaX10'] = x10 - x9
            
            # DeltaY10: normalizzato
            if vdist > 1e-6:
                p['DeltaY10'] = (y10 - y9) / vdist
            
            # DeltaX11
            p['DeltaX11'] = abs(x12 - x11)

            if self.slat_type.get() != "None":
                # Slat parameters
                x13, y13 = self.main_slot_slat_control_points[0]
                x14, y14 = self.main_slot_slat_control_points[1]
                x15, y15 = self.main_slot_slat_control_points[2]
                x16, y16 = self.main_slot_slat_control_points[3]
                x17, y17 = self.main_slot_slat_control_points[4]
                x18, y18 = self.main_slot_slat_control_points[5]
                x19, y19 = self.slat_control_points[0]
                x20, y20 = self.slat_control_points[1]
                x21, y21 = self.slat_control_points[2]
                x22, y22 = self.slat_control_points[3]
                x23, y23 = self.slat_control_points[4]
                x24, y24 = self.slat_control_points[5]

                # Direct coordinate updates
                p['x13'] = x13
                p['y13'] = y13
                p['x14'] = x14
                p['y14'] = y14
                p['x15'] = x15
                p['y15'] = y15
                p['x16'] = x16
                p['y16'] = y16
                p['x17'] = x17
                p['y17'] = y17
                p['x18'] = x18
                p['y18'] = y18
                p['x19'] = x19
                p['y19'] = y19
                p['x20'] = x20
                p['y20'] = y20
                p['x21'] = x21
                p['y21'] = y21
                p['x22'] = x22
                p['y22'] = y22
                p['x23'] = x23
                p['y23'] = y23
                p['x24'] = x24
                p['y24'] = y24
                
                # Hinge point for slat
                p['xh_slat'] = self.slat_hinge_point[0]
                p['yh_slat'] = self.slat_hinge_point[1]
            
            # Update parameter GUI fields
            # for param, var in self.param_vars.items():
            #     if param in p:
            #         var.set(round(p[param], 6))

            self.update_tree_parameters()
                    
        except Exception as e:
            self.status_label.config(text=f"Error updating parameters: {str(e)}")

    def restart_code(self):
        """Restart the application"""
        self.root.destroy()
        self.__init__(tk.Tk())

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

                if self.flap_type.get() != "None" and self.slat_type.get() != "None":
                    self.generate_flap()
                    self.generate_slat()
                    self.assemble_components()
                    self.plot_results()
                    self.status_label.config(text="Flap and Slat generated! Drag control points to modify")
                    self.calculate_gap_overlap()

                elif self.flap_type.get() != "None" and self.slat_type.get() == "None":
                    self.generate_flap()
                    self.assemble_components()
                    self.plot_results()
                    self.status_label.config(text="Flap generated! Drag control points to modify")
                    self.calculate_gap_overlap()
                
                elif self.flap_type.get() == "None" and self.slat_type.get() != "None":
                    self.generate_slat()
                    self.assemble_components()
                    self.plot_results()
                    self.status_label.config(text="Slat generated! Drag control points to modify")
                    self.calculate_gap_overlap()
                
                elif self.flap_type.get() == "None" and self.slat_type.get() == "None":
                    self.assemble_components()
                    self.plot_results()
                    self.status_label.config(text="Basic airfoil loaded successfully!")

                

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

        self.f_upper = interp1d(self.x_upper, self.y_upper, kind='linear', bounds_error=False, fill_value='extrapolate')
        self.f_lower = interp1d(self.x_lower, self.y_lower, kind='linear', bounds_error=False, fill_value='extrapolate')
    
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
    
    def generate_slat(self):
        """Generate slat based on current parameters"""
        if self.airfoil_data is None:
            self.status_label.config(text="Please load an airfoil first!")
            return
        
        try:
            self.generate_main_slot_slat_control_points()
            self.generate_slat_control_points()
            self.slat_hinge_point = np.array([self.params['xh_slat'], self.params['yh_slat']])
            self.main_slot_slat_curve = self.bezier_curve(self.main_slot_slat_control_points)
            self.slat_curve = self.bezier_curve(self.slat_control_points)

        except Exception as e:
            self.status_label.config(text=f"Error generating slat: {str(e)}")

    def generate_main_slot_slat_control_points(self):
        """Generate main slot slat control points"""
        p = self.params

        x13, y13 = p['x13'], p['y13']
        x14, y14 = p['x14'], p['y14']
        x15, y15 = p['x15'], p['y15']
        x16, y16 = p['x16'], p['y16']
        x17, y17 = p['x17'], p['y17']
        x18, y18 = p['x18'], p['y18']
        
        self.main_slot_slat_control_points = np.array([[x13, y13], [x14, y14], [x15, y15], [x16, y16], [x17, y17], [x18, y18]])
    
    def generate_slat_control_points(self):
        """Generate slat control points"""
        p = self.params

        x19, y19 = p['x19'], p['y19']
        x20, y20 = p['x20'], p['y20']
        x21, y21 = p['x21'], p['y21']
        x22, y22 = p['x22'], p['y22']
        x23, y23 = p['x23'], p['y23']
        x24, y24 = p['x24'], p['y24']

        self.slat_control_points = np.array([[x19, y19], [x20, y20], [x21, y21], [x22, y22], [x23, y23], [x24, y24]])

    def generate_flap(self):
        """Generate flap based on current parameters"""
        if self.airfoil_data is None:
            self.status_label.config(text="Please load an airfoil first!")
            return
        
        try:
            self.generate_main_slot_flap_control_points()
            self.generate_flap_control_points()
            self.flap_hinge_point = np.array([self.params['xh_flap'], self.params['yh_flap']])
            self.main_slot_flap_curve = self.bezier_curve(self.main_slot_flap_control_points)
            self.flap_curve = self.bezier_curve(self.flap_control_points)
        except Exception as e:
            self.status_label.config(text=f"Error generating flap: {str(e)}")
        
    def generate_main_slot_flap_control_points(self):
        """Generate main slot flap control points"""
        p = self.params

        # Point calculations (same as original)
        x1, y1 = p['x1'], self.f_lower(p['x1'])
        x2p = x1 - p['Delta_X1']
        y2p = self.f_lower(x2p)
        slope_lower = np.arctan((y2p - y1) / (x2p - x1))
        x2 = p['x2']
        y2 = y1 + slope_lower * (x2 - x1)
        
        x3 = p['x3']
        y3p = self.f_upper(x3)
        delta_y3 = p['FractionY'] * (y3p - y2)
        y3 = y2 + delta_y3
        
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
        
        self.main_slot_flap_control_points = np.array([[x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5], [x6, y6]])
    
    def generate_flap_control_points(self):
        """Generate flap control points"""
        p = self.params

        x7, y7 = p['x7'], self.f_upper(p['x7'])
        x8, y8 = p['x6'], self.f_upper(p['x6'])
        
        x9 = self.main_slot_flap_control_points[2, 0] + p['DeltaX9']
        vdist_9 = abs(self.f_upper(p['x6']) - self.f_lower(p['x6']))
        y9 = self.main_slot_flap_control_points[2, 1] + p['DeltaY9'] * vdist_9
        
        x10 = x9 + p['DeltaX10']
        y10 = y9 + p['DeltaY10'] * vdist_9
        
        x12, y12 = p['x12'], self.f_lower(p['x12'])
        x11 = x12 - p['DeltaX11']
        x12_bis = x12 + 0.00001
        y12_bis = self.f_lower(x12_bis)
        slope_y11 = np.arctan((y12_bis - y12) / (x12_bis - x12))
        y11 = y12 - slope_y11 * (x12 - x11)
        
        self.flap_control_points = np.array([[x7, y7], [x8, y8], [x9, y9], [x10, y10], [x11, y11], [x12, y12]])
    
    def assemble_components(self):
        """Assemble main and flap components"""

        if self.airfoil_data is None:
            self.status_label.config(text="Please load an airfoil first!")
            return

        if self.flap_type.get() != "None" and self.slat_type.get() == "None":
            p = self.params
            
            ix6 = np.argmin(np.abs(p['x6'] - self.x_upper))
            ix1 = np.argmin(np.abs(p['x1'] - self.x_lower))
            ix7 = np.argmin(np.abs(p['x7'] - self.x_upper))
            ix12 = np.argmin(np.abs(p['x12'] - self.x_lower))
            
            # Main component
            main_x = np.concatenate([self.x_upper[ix6:], self.x_lower[:ix1], self.main_slot_flap_curve[:, 0]])
            main_y = np.concatenate([self.y_upper[ix6:], self.y_lower[:ix1], self.main_slot_flap_curve[:, 1]])
            
            cos_twist = np.cos(-np.deg2rad(p['Twist Main']))
            sin_twist = np.sin(-np.deg2rad(p['Twist Main']))
            self.main_component_undeflected = np.column_stack([main_x, main_y])
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
            self.flap_component_undeflected = np.column_stack([flap_x, flap_y])
            self.slat_component = None  # No slat component

        elif self.flap_type.get() != "None" and self.slat_type.get() != "None":
            p = self.params
            
            ix6 = np.argmin(np.abs(p['x6'] - self.x_upper))
            ix1 = np.argmin(np.abs(p['x1'] - self.x_lower))
            ix7 = np.argmin(np.abs(p['x7'] - self.x_upper))
            ix12 = np.argmin(np.abs(p['x12'] - self.x_lower))
            ix13 = np.argmin(np.abs(p['x13'] - self.x_upper))
            ix18 = np.argmin(np.abs(p['x18'] - self.x_lower))
            ix19 = np.argmin(np.abs(p['x19'] - self.x_lower))
            ix24 = np.argmin(np.abs(p['x24'] - self.x_upper))
            
            # Main component
            main_x = np.concatenate([self.x_upper[ix6:ix13], self.main_slot_slat_curve[:, 0], self.x_lower[ix18+1:ix1], self.main_slot_flap_curve[:, 0]])
            main_y = np.concatenate([self.y_upper[ix6:ix13], self.main_slot_slat_curve[:, 1], self.y_lower[ix18+1:ix1], self.main_slot_flap_curve[:, 1]])
            
            cos_twist = np.cos(-np.deg2rad(p['Twist Main']))
            sin_twist = np.sin(-np.deg2rad(p['Twist Main']))
            self.main_component_undeflected = np.column_stack([main_x, main_y])
            self.main_component = np.column_stack([
                main_x * cos_twist - main_y * sin_twist,
                main_x * sin_twist + main_y * cos_twist
            ])
            
            # Flap component (undeflected)
            flap_x = np.concatenate([self.x_upper[:ix7], self.flap_curve[:, 0], self.x_lower[ix12+1:]])
            flap_y = np.concatenate([self.y_upper[:ix7], self.flap_curve[:, 1], self.y_lower[ix12+1:]])
            
            # Apply main rotation only for display
            self.flap_component_undeflected = np.column_stack([flap_x, flap_y])
            self.flap_component = np.column_stack([
                flap_x * cos_twist - flap_y * sin_twist,
                flap_x * sin_twist + flap_y * cos_twist
            ])

            # Slat component
            slat_x = np.concatenate([self.x_upper[ix24:-1], self.x_lower[:ix19], self.slat_curve[:, 0]])
            slat_y = np.concatenate([self.y_upper[ix24:-1], self.y_lower[:ix19], self.slat_curve[:, 1]])
            slat_x[0] = p['x24']  
            slat_y[0] = self.f_upper(p['x24']) 

            # Apply main rotation only for display
            self.slat_component_undeflected = np.column_stack([slat_x, slat_y])
            self.slat_component = np.column_stack([
                slat_x * cos_twist - slat_y * sin_twist,
                slat_x * sin_twist + slat_y * cos_twist
            ])

        elif self.flap_type.get() == "None" and self.slat_type.get() != "None":

            p = self.params
            
            ix13 = np.argmin(np.abs(p['x13'] - self.x_upper))
            ix18 = np.argmin(np.abs(p['x18'] - self.x_lower))
            ix19 = np.argmin(np.abs(p['x19'] - self.x_lower))
            ix24 = np.argmin(np.abs(p['x24'] - self.x_upper))
            
            # Main component
            main_x = np.concatenate([self.x_upper[:ix13], self.main_slot_slat_curve[:, 0], self.x_lower[ix18+1:]])
            main_y = np.concatenate([self.y_upper[:ix13], self.main_slot_slat_curve[:, 1], self.y_lower[ix18+1:]])
            
            cos_twist = np.cos(-np.deg2rad(p['Twist Main']))
            sin_twist = np.sin(-np.deg2rad(p['Twist Main']))
            self.main_component_undeflected = np.column_stack([main_x, main_y])
            self.main_component = np.column_stack([
                main_x * cos_twist - main_y * sin_twist,
                main_x * sin_twist + main_y * cos_twist
            ])
            
            # Slat component
            slat_x = np.concatenate([self.x_upper[ix24:-1], self.x_lower[:ix19], self.slat_curve[:, 0]])
            slat_y = np.concatenate([self.y_upper[ix24:-1], self.y_lower[:ix19], self.slat_curve[:, 1]])
            slat_x[0] = p['x24']  
            slat_y[0] = self.f_upper(p['x24']) 

            # Apply main rotation only for display
            self.slat_component_undeflected = np.column_stack([slat_x, slat_y])
            self.slat_component = np.column_stack([
                slat_x * cos_twist - slat_y * sin_twist,
                slat_x * sin_twist + slat_y * cos_twist
            ])
            self.flap_component = None

        elif self.flap_type.get() == "None" and self.slat_type.get() == "None":
            # No flap or slat, just the main component
            p = self.params

            # Main component
            main_x = np.concatenate([self.x_upper, self.x_lower])
            main_y = np.concatenate([self.y_upper, self.y_lower])
            
            cos_twist = np.cos(-np.deg2rad(p['Twist Main']))
            sin_twist = np.sin(-np.deg2rad(p['Twist Main']))
            self.main_component_undeflected = np.column_stack([main_x, main_y])
            self.main_component = np.column_stack([
                main_x * cos_twist - main_y * sin_twist,
                main_x * sin_twist + main_y * cos_twist
            ])
            
            self.flap_component = None
            self.slat_component = None

    
    def get_deflected_flap(self):
        """Get flap component with deflection applied"""
        p = self.params
        
        # Apply deflection
        xh, yh = p['xh_flap'], p['yh_flap']
        delta = np.deg2rad(p['Flap Deflection'])
        
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
        """Plot results flap"""
        self.ax.clear()

        if self.airfoil_data is None:
            self.status_label.config(text="Please load an airfoil first!")
            return
        
        self.ax.plot(self.airfoil_data[:, 0], self.airfoil_data[:, 1], 'k-', linewidth=1, alpha=0.5)

        if self.params['Twist Main'] != 0:
            self.ax.plot(self.main_component_undeflected[:, 0], self.main_component_undeflected[:, 1], 'k-', linewidth=1, alpha=0.5, label='Original')
            if self.flap_component_undeflected is not None:
                self.ax.plot(self.flap_component_undeflected[:, 0], self.flap_component_undeflected[:, 1], 'k-', linewidth=1, alpha=0.5)
            if self.slat_component_undeflected is not None:
                self.ax.plot(self.slat_component_undeflected[:, 0], self.slat_component_undeflected[:, 1], 'k-', linewidth=1, alpha=0.5)
        
        if self.main_component is not None:
            self.main_line = self.ax.plot(self.main_component[:, 0], self.main_component[:, 1], 'b-', linewidth=2.5, label='Main')[0]
        
        if self.flap_component is not None:
            # Plot deflected flap if deflection is not zero
            if abs(self.params['Flap Deflection']) > 1e-6:
                deflected_flap = self.apply_flap_deflection()
                self.flap_line = self.ax.plot(self.flap_component[:, 0], self.flap_component[:, 1], 'r--', linewidth=1.5, alpha=0.5, label='Flap (undeflected)')[0]
                self.deflected_flap_line = self.ax.plot(deflected_flap[:, 0], deflected_flap[:, 1], 'r-', linewidth=2.5, label=f'Flap (deflected {self.params["Flap Deflection"]:.1f}°)')[0]
            else:
                self.flap_line = self.ax.plot(self.flap_component[:, 0], self.flap_component[:, 1], 'r-', linewidth=1.5, label='Flap (undeflected)')[0]

        if self.slat_component is not None:
            # Ensures main il plotted once
            if self.flap_component is None:
                self.main_slot_slat_curve_line = self.ax.plot(self.main_slot_slat_curve[:, 0], self.main_slot_slat_curve[:, 1], 'b:', linewidth=1.5, alpha=0.7, label='Main Curve')[0]
        
            # Plot deflected slat if deflection is not zero
            if abs(self.params['Slat Deflection']) > 1e-6:
                deflected_slat = self.apply_slat_deflection()
                self.slat_line = self.ax.plot(self.slat_component[:, 0], self.slat_component[:, 1], 'g--', linewidth=1.5, alpha=0.5, label='Slat (undeflected)')[0]
                self.deflected_slat_line = self.ax.plot(deflected_slat[:, 0], deflected_slat[:, 1], 'g-', linewidth=2.5, label=f'Slat (deflected {self.params["Slat Deflection"]:.1f}°)')[0]
            else:
                self.slat_line = self.ax.plot(self.slat_component[:, 0], self.slat_component[:, 1], 'g-', linewidth=1.5, label='Slat (undeflected)')[0]

        self.create_interactive_points()
        
        self.ax.set_xlabel('x/c')
        self.ax.set_ylabel('y/c')
        self.ax.grid(True, alpha=0.3)
        self.ax.axis('equal')
        self.ax.legend()
        self.ax.set_title('Interactive Airfoil-Flap-Slat Designer')
        
        self.canvas.draw()
    
    def create_interactive_points(self):
        """Create interactive control points with better hit detection"""
        # Clear existing artists and annotations
        for artist in self.main_slot_flap_artists + self.flap_artists + self.main_slot_slat_artists + self.slat_artists:
            if artist in self.ax.collections:
                artist.remove()
        if self.flap_hinge_artist and self.flap_hinge_artist in self.ax.collections:
            self.flap_hinge_artist.remove()
        if self.slat_hinge_artist and self.slat_hinge_artist in self.ax.collections:
            self.slat_hinge_artist.remove()

        # Clear existing control lines
        if hasattr(self, 'main_slot_flap_control_lines'):
            for line in self.main_slot_flap_control_lines:
                line.remove()
        if hasattr(self, 'flap_control_lines'):
            for line in self.flap_control_lines:
                line.remove()
        if hasattr(self, 'main_slot_slat_control_lines'):
            for line in self.main_slot_slat_control_lines:
                line.remove()
        if hasattr(self, 'slat_control_lines'):
            for line in self.slat_control_lines:
                line.remove()

        # RIMUOVERE le annotazioni esistenti
        for ann in self.point_annotations:
            ann.remove()

        self.main_slot_flap_artists = []
        self.flap_artists = []
        self.main_slot_flap_control_lines = []
        self.flap_control_lines = []
        self.point_annotations = []  
        self.main_slot_slat_artists = []
        self.slat_artists = []
        self.main_slot_slat_control_lines = []
        self.slat_control_lines = []

        if self.flap_type.get() != "None":
            # Create slot control lines
            if self.main_slot_flap_control_points is not None:
                for i in range(len(self.main_slot_flap_control_points) - 1):
                    x_data = [self.main_slot_flap_control_points[i][0], self.main_slot_flap_control_points[i+1][0]]
                    y_data = [self.main_slot_flap_control_points[i][1], self.main_slot_flap_control_points[i+1][1]]
                    line = self.ax.plot(x_data, y_data, 'b--', linewidth=0.8, alpha=0.6, zorder=5)[0]
                    self.main_slot_flap_control_lines.append(line)
            
            # Create flap control lines
            if self.flap_control_points is not None:
                for i in range(len(self.flap_control_points) - 1):
                    x_data = [self.flap_control_points[i][0], self.flap_control_points[i+1][0]]
                    y_data = [self.flap_control_points[i][1], self.flap_control_points[i+1][1]]
                    line = self.ax.plot(x_data, y_data, 'r--', linewidth=0.8, alpha=0.6, zorder=5)[0]
                    self.flap_control_lines.append(line)
            
            # Create slot control points
            if self.main_slot_flap_control_points is not None:
                for i, point in enumerate(self.main_slot_flap_control_points):
                    artist = self.ax.scatter(point[0], point[1], c='blue', s=80, marker='o', 
                                        edgecolor='darkblue', linewidth=2, zorder=15, alpha=0.9)
                    self.main_slot_flap_artists.append(artist)
            
            # Create flap control points
            if self.flap_control_points is not None:
                for i, point in enumerate(self.flap_control_points):
                    artist = self.ax.scatter(point[0], point[1], c='red', s=80, marker='s', 
                                        edgecolor='darkred', linewidth=2, zorder=15, alpha=0.9)
                    self.flap_artists.append(artist)
            
            # Create hinge point
            if self.flap_hinge_point is not None:
                self.flap_hinge_artist = self.ax.scatter(self.flap_hinge_point[0], self.flap_hinge_point[1], 
                                                c='green', s=100, marker='D', 
                                                edgecolor='darkgreen', linewidth=2, zorder=15, alpha=0.9)

        if self.slat_type.get() != "None":
            # Create slot control lines for slat
            if self.main_slot_slat_control_points is not None:
                for i in range(len(self.main_slot_slat_control_points) - 1):
                    x_data = [self.main_slot_slat_control_points[i][0], self.main_slot_slat_control_points[i+1][0]]
                    y_data = [self.main_slot_slat_control_points[i][1], self.main_slot_slat_control_points[i+1][1]]
                    line = self.ax.plot(x_data, y_data, 'b--', linewidth=0.8, alpha=0.6, zorder=5)[0]
                    self.main_slot_slat_control_lines.append(line)

            # Create slat control lines
            if self.slat_control_points is not None:
                for i in range(len(self.slat_control_points) - 1):
                    x_data = [self.slat_control_points[i][0], self.slat_control_points[i+1][0]]
                    y_data = [self.slat_control_points[i][1], self.slat_control_points[i+1][1]]
                    line = self.ax.plot(x_data, y_data, 'r--', linewidth=0.8, alpha=0.6, zorder=5)[0]
                    self.slat_control_lines.append(line)

            # Create slot control points
            if self.main_slot_slat_control_points is not None:
                for i, point in enumerate(self.main_slot_slat_control_points):
                    artist = self.ax.scatter(point[0], point[1], c='blue', s=80, marker='o', 
                                        edgecolor='darkblue', linewidth=2, zorder=15, alpha=0.9)
                    self.main_slot_slat_artists.append(artist)
            
            # Create slat control points
            if self.slat_control_points is not None:
                for i, point in enumerate(self.slat_control_points):
                    artist = self.ax.scatter(point[0], point[1], c='green', s=80, marker='s', 
                                        edgecolor='darkgreen', linewidth=2, zorder=15, alpha=0.9)
                    self.slat_artists.append(artist)
            
            # Create hinge point
            if self.slat_hinge_point is not None:
                self.slat_hinge_artist = self.ax.scatter(self.slat_hinge_point[0], self.slat_hinge_point[1], 
                                                c='orange', s=100, marker='D', 
                                                edgecolor='darkorange', linewidth=2, zorder=15, alpha=0.9)




    def update_curves_and_components(self):
        """Update curves and components"""
        try:
            if self.flap_component is not None:
                self.main_slot_flap_curve = self.bezier_curve(self.main_slot_flap_control_points)
                self.flap_curve = self.bezier_curve(self.flap_control_points)

            if self.slat_component is not None:
                self.main_slot_slat_curve = self.bezier_curve(self.main_slot_slat_control_points)
                self.slat_curve = self.bezier_curve(self.slat_control_points)
            
            # Update parameter values
            self.params['xh_flap'] = self.flap_hinge_point[0]
            self.params['yh_flap'] = self.flap_hinge_point[1]
            self.param_vars['xh_flap'].set(self.flap_hinge_point[0])
            self.param_vars['yh_flap'].set(self.flap_hinge_point[1])

            self.params['xh_slat'] = self.slat_hinge_point[0]
            self.params['yh_slat'] = self.slat_hinge_point[1]
            self.param_vars['xh_slat'].set(self.slat_hinge_point[0])
            self.param_vars['yh_slat'].set(self.slat_hinge_point[1])
            
            self.assemble_components()
            
            # Update plot lines
            if self.main_line:
                self.main_line.set_data(self.main_component[:, 0], self.main_component[:, 1])
            if self.flap_line:
                self.flap_line.set_data(self.flap_component[:, 0], self.flap_component[:, 1])
            if self.main_slot_flap_curve_line:
                self.main_slot_flap_curve_line.set_data(self.main_slot_flap_curve[:, 0], self.main_slot_flap_curve[:, 1])

            if self.flap_curve_line:
                self.flap_curve_line.set_data(self.flap_curve[:, 0], self.flap_curve[:, 1])
            if self.slat_line:
                self.slat_line.set_data(self.slat_component[:, 0], self.slat_component[:, 1])
            
            if self.main_slot_slat_curve_line:
                self.main_slot_slat_curve_line.set_data(self.main_slot_slat_curve[:, 0], self.main_slot_slat_curve[:, 1])
            
            # Update control lines
            if hasattr(self, 'main_slot_flap_curve_line') and self.main_slot_flap_control_points is not None:
                for i, line in enumerate(self.main_slot_flap_control_lines):
                    if i < len(self.main_slot_flap_control_points) - 1:
                        x_data = [self.main_slot_flap_control_points[i][0], self.main_slot_flap_control_points[i+1][0]]
                        y_data = [self.main_slot_flap_control_points[i][1], self.main_slot_flap_control_points[i+1][1]]
                        line.set_data(x_data, y_data)

            if hasattr(self, 'main_slot_slat_curve_line') and self.main_slot_slat_control_points is not None:
                for i, line in enumerate(self.main_slot_slat_control_lines):
                    if i < len(self.main_slot_slat_control_points) - 1:
                        x_data = [self.main_slot_slat_control_points[i][0], self.main_slot_slat_control_points[i+1][0]]
                        y_data = [self.main_slot_slat_control_points[i][1], self.main_slot_slat_control_points[i+1][1]]
                        line.set_data(x_data, y_data)
            
            if hasattr(self, 'flap_control_lines') and self.flap_control_points is not None:
                for i, line in enumerate(self.flap_control_lines):
                    if i < len(self.flap_control_points) - 1:
                        x_data = [self.flap_control_points[i][0], self.flap_control_points[i+1][0]]
                        y_data = [self.flap_control_points[i][1], self.flap_control_points[i+1][1]]
                        line.set_data(x_data, y_data)

            if hasattr(self, 'slat_control_lines') and self.slat_control_points is not None:
                for i, line in enumerate(self.slat_control_lines):
                    if i < len(self.slat_control_points) - 1:
                        x_data = [self.slat_control_points[i][0], self.slat_control_points[i+1][0]]
                        y_data = [self.slat_control_points[i][1], self.slat_control_points[i+1][1]]
                        line.set_data(x_data, y_data)
            
            # Update deflected flap line if it exists
            if hasattr(self, 'deflected_flap_line') and abs(self.params['Flap Deflection']) > 1e-6:
                deflected_flap = self.apply_flap_deflection()
                self.deflected_flap_line.set_data(deflected_flap[:, 0], deflected_flap[:, 1])
                # Update label
                self.deflected_flap_line.set_label(f'Flap (deflected {self.params["Flap Deflection"]:.1f}°)')
                self.calculate_gap_overlap()
            elif hasattr(self, 'deflected_flap_line') and abs(self.params['Flap Deflection']) <= 1e-6:
                # Remove deflected flap line if deflection is zero
                self.deflected_flap_line.remove()
                delattr(self, 'deflected_flap_line')
                self.calculate_gap_overlap()
            elif not hasattr(self, 'deflected_flap_line') and abs(self.params['Flap Deflection']) > 1e-6:
                # Create deflected flap line if it doesn't exist
                deflected_flap = self.apply_flap_deflection()
                self.deflected_flap_line = self.ax.plot(deflected_flap[:, 0], deflected_flap[:, 1], 'r-', linewidth=2.5, 
                                                    label=f'Flap (deflected {self.params["Flap Deflection"]:.1f}°)')[0]
                self.calculate_gap_overlap()
                
            if hasattr(self, 'deflected_slat_line') and abs(self.params['Slat Deflection']) > 1e-6:
                deflected_slat = self.apply_slat_deflection()
                self.deflected_slat_line.set_data(deflected_slat[:, 0], deflected_slat[:, 1])
                # Update label
                self.deflected_slat_line.set_label(f'Slat (deflected {self.params["Slat Deflection"]:.1f}°)')
            elif hasattr(self, 'deflected_slat_line') and abs(self.params['Slat Deflection']) <= 1e-6:
                # Remove deflected slat line if deflection is zero
                self.deflected_slat_line.remove()
                delattr(self, 'deflected_slat_line')
            elif not hasattr(self, 'deflected_slat_line') and abs(self.params['Slat Deflection']) > 1e-6:
                # Create deflected slat line if it doesn't exist
                deflected_slat = self.apply_slat_deflection()
                self.deflected_slat_line = self.ax.plot(deflected_slat[:, 0], deflected_slat[:, 1], 'g-', linewidth=2.5, 
                                                    label=f'Slat (deflected {self.params["Slat Deflection"]:.1f}°)')[0]
            
            # Update point annotations
            for i, ann in enumerate(self.point_annotations):
                if i < len(self.main_slot_flap_control_points):
                    ann.set_position(self.main_slot_flap_control_points[i])
                elif i < len(self.main_slot_flap_control_points) + len(self.flap_control_points):
                    idx = i - len(self.main_slot_flap_control_points)
                    ann.set_position(self.flap_control_points[idx])
                else:
                    ann.set_position(self.flap_hinge_point)

                if i < len(self.main_slot_slat_control_points):
                    ann.set_position(self.main_slot_slat_control_points[i])
                elif i < len(self.main_slot_slat_control_points) + len(self.slat_control_points):
                    idx = i - len(self.main_slot_slat_control_points)
                    ann.set_position(self.slat_control_points[idx])
                else:
                    ann.set_position(self.slat_hinge_point)

            
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
            
            # Check main slot flap control points
            if self.main_slot_flap_control_points is not None:
                for i, point in enumerate(self.main_slot_flap_control_points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance:
                        candidates.append((distance, 'main_slot_flap', i, point))
            
            # Check flap control points
            if self.flap_control_points is not None:
                for i, point in enumerate(self.flap_control_points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance:
                        candidates.append((distance, 'flap', i, point))
            
            # Check hinge point
            if self.flap_hinge_point is not None:
                distance = np.sqrt((event.xdata - self.flap_hinge_point[0])**2 + (event.ydata - self.flap_hinge_point[1])**2)
                if distance < tolerance:
                    candidates.append((distance, 'hinge_flap', 0, self.flap_hinge_point))

            # Check slat control points 
            if self.main_slot_slat_control_points is not None:
                for i, point in enumerate(self.main_slot_slat_control_points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance:
                        candidates.append((distance, 'main_slot_slat', i, point))

            # Check slat flap control points
            if self.slat_control_points is not None:
                for i, point in enumerate(self.slat_control_points):
                    distance = np.sqrt((event.xdata - point[0])**2 + (event.ydata - point[1])**2)
                    if distance < tolerance:
                        candidates.append((distance, 'slat', i, point))

            # Check hinge point for slat
            if self.slat_hinge_point is not None:
                distance = np.sqrt((event.xdata - self.slat_hinge_point[0])**2 + (event.ydata - self.slat_hinge_point[1])**2)
                if distance < tolerance:
                    candidates.append((distance, 'hinge_slat', 0, self.slat_hinge_point))
            
            # Selezionare il punto più vicino
            if candidates:
                # Ordinare per distanza (più vicino prima)
                candidates.sort(key=lambda x: x[0])
                distance, point_type, point_idx, point = candidates[0]
                
                self.dragging_point = (point_type, point_idx)
                self.drag_offset = (event.xdata - point[0], event.ydata - point[1])
                
                if point_type == 'main_slot_flap':
                    self.status_label.config(text=f"Dragging Main Point {point_idx+1}")
                elif point_type == 'flap':
                    self.status_label.config(text=f"Dragging Flap Point {point_idx+1}")
                elif point_type == 'hinge_flap':
                    self.status_label.config(text="Dragging Hinge Point")
                elif point_type == 'main_slot_slat':
                    self.status_label.config(text=f"Dragging Main Slot Slat Point {point_idx+1}")
                elif point_type == 'slat':
                    self.status_label.config(text=f"Dragging Slat Point {point_idx+1}")
                elif point_type == 'hinge_slat':
                    self.status_label.config(text="Dragging Hinge Slat Point")
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
            if curve_type == 'main_slot_flap':
                self.main_slot_flap_control_points[point_idx] = [new_x, new_y]
                # Update corresponding scatter plot
                self.main_slot_flap_artists[point_idx].set_offsets([[new_x, new_y]])
            elif curve_type == 'flap':
                self.flap_control_points[point_idx] = [new_x, new_y]
                # Update corresponding scatter plot
                self.flap_artists[point_idx].set_offsets([[new_x, new_y]])
            elif curve_type == 'hinge_flap':
                self.flap_hinge_point = np.array([new_x, new_y])
                # Update corresponding scatter plot
                self.flap_hinge_artist.set_offsets([[new_x, new_y]])
            elif curve_type == 'main_slot_slat':
                self.main_slot_slat_control_points[point_idx] = [new_x, new_y]
                # Update corresponding scatter plot
                self.main_slot_slat_artists[point_idx].set_offsets([[new_x, new_y]])
            elif curve_type == 'slat':
                self.slat_control_points[point_idx] = [new_x, new_y]
                # Update corresponding scatter plot
                self.slat_artists[point_idx].set_offsets([[new_x, new_y]])
            elif curve_type == 'hinge_slat':
                self.slat_hinge_point = np.array([new_x, new_y])
                # Update corresponding scatter plot
                self.slat_hinge_artist.set_offsets([[new_x, new_y]])
            
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
                    'main_slot_flap_control_points': self.main_slot_flap_control_points.tolist() if self.main_slot_flap_control_points is not None else None,
                    'flap_control_points': self.flap_control_points.tolist() if self.flap_control_points is not None else None,
                    'hinge_point_flap': self.flap_hinge_point.tolist() if self.flap_hinge_point is not None else None,
                    'main_slot_slat_control_points': self.main_slot_slat_control_points.tolist() if self.main_slot_slat_control_points is not None else None,
                    'slat_control_points': self.slat_control_points.tolist() if self.slat_control_points is not None else None,
                    'hinge_point_slat': self.slat_hinge_point.tolist() if self.slat_hinge_point is not None else None,
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
                if config.get('main_slot_flap_control_points') is not None:
                    self.main_slot_flap_control_points = np.array(config['main_slot_flap_control_points'])
                if config.get('flap_control_points') is not None:
                    self.flap_control_points = np.array(config['flap_control_points'])
                if config.get('hinge_point_flap') is not None:
                    self.flap_hinge_point = np.array(config['hinge_point_flap'])
                if config.get('main_slot_slat_control_points') is not None:
                    self.main_slot_slat_control_points = np.array(config['main_slot_slat_control_points'])
                if config.get('slat_control_points') is not None:
                    self.slat_control_points = np.array(config['slat_control_points'])
                if config.get('hinge_point_slat') is not None:
                    self.slat_hinge_point = np.array(config['hinge_point_slat'])

                if 'magnetic_snap' in config:
                    self.magnetic_snap = config['magnetic_snap']
                    if self.magnetic_snap:
                        self.magnet_style.configure("Active.TButton", background="lightgreen")
                        self.magnet_button.config(style="Active.TButton")
                    else:
                        self.magnet_button.config(style="TButton")
                
                # Update curves first
                self.main_slot_flap_curve = self.bezier_curve(self.main_slot_flap_control_points)
                self.flap_curve = self.bezier_curve(self.flap_control_points)
                self.main_slot_slat_curve = self.bezier_curve(self.main_slot_slat_control_points)
                self.slat_curve = self.bezier_curve(self.slat_control_points)
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