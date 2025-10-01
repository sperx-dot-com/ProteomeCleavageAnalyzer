import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import re 
from sequence_analyzer import SequenceAnalyzer, GOFilter
from structure_analyzer import StructureAnalyzer
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.styles import Font
import webbrowser

class SequenceAnalyzerUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Sequence Analyzer")
        self.root.geometry("1500x1000")

        # Global paths
        self.GENOME_FOLDER_PATH = "genomes"
        self.GO_OBO_PATH = os.path.join("go_files", "go.obo")
        
        # Initialize business logic classes
        self.sequence_analyzer = SequenceAnalyzer()
        self.structure_analyzer = StructureAnalyzer()
        self.go_filter = GOFilter()
        
        # Make window resizable
        self.root.rowconfigure(4, weight=1)  # Update row number for results frame
        self.root.columnconfigure(0, weight=1)
        
        # Current hits (global state)
        self.current_hits = []
        self.filtered_hits = []  # Add this new attribute to store filtered hits
        # Dictionary to map tree item ids to hit indexes
        self.tree_id_to_hit_index = {}

        self.current_go_filter = ""

        # Initialize UI
        self.setup_top_frame()
        self.setup_param_frame()
        self.setup_search_button()
        self.setup_go_filter_frame()  # Move this below search button
        self.setup_results_frame()
        self.setup_bottom_frame()

    def setup_top_frame(self):
        # ----- Top Frame: Genome File Selection ----- #
        top_frame = tk.Frame(self.root)
        top_frame.grid(row=0, column=0, pady=10, padx=10, sticky="ew")
        
        tk.Label(top_frame, text="Genome File:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.fasta_file_entry = tk.Entry(top_frame, width=50)
        self.fasta_file_entry.grid(row=0, column=1, padx=5)
        
        default_genomes = self.get_default_genomes()
        self.genome_var = tk.StringVar()
        self.genome_dropdown = ttk.Combobox(top_frame, textvariable=self.genome_var, 
                                           values=default_genomes, state="readonly", width=20)
        self.genome_dropdown.grid(row=0, column=2, padx=5)
        self.genome_dropdown.bind("<<ComboboxSelected>>", self.update_fasta_entry)
        
        tk.Button(top_frame, text="Browse", 
                 command=lambda: self.browse_file(self.fasta_file_entry, 
                                                [("FASTA Files", "*.fasta *.fa *.faa *.txt")])).grid(row=0, column=3, padx=5)

    def setup_param_frame(self):
        # ----- Middle Frame: Parameters ----- #
        param_frame = tk.Frame(self.root)
        param_frame.grid(row=1, column=0, pady=10, padx=10, sticky="ew")
        
        tk.Label(param_frame, text="Regex Pattern:").grid(row=0, column=0, sticky=tk.W, padx=5)
        self.regex_entry = tk.Entry(param_frame, width=30)
        self.regex_entry.grid(row=0, column=1, padx=5)
        
        tk.Label(param_frame, text="P1' Pattern (comma separated):").grid(row=0, column=2, sticky=tk.W, padx=5)
        self.p1_entry = tk.Entry(param_frame, width=15)
        self.p1_entry.grid(row=0, column=3, padx=5)
        
        self.at_cterm_var = tk.BooleanVar()
        tk.Checkbutton(param_frame, text="At C-Terminus", 
                      variable=self.at_cterm_var).grid(row=1, column=1, sticky=tk.W, padx=5)
        self.close_cterm_var = tk.BooleanVar()
        tk.Checkbutton(param_frame, text="Close to C-Terminus (within 10 AA)", 
                      variable=self.close_cterm_var).grid(row=1, column=2, sticky=tk.W, padx=5)
        
        tk.Label(param_frame, text="GO Annotation (GAF):").grid(row=2, column=0, sticky=tk.W, padx=5, pady=(10,0))
        self.go_entry = tk.Entry(param_frame, width=30)
        self.go_entry.grid(row=2, column=1, padx=5, pady=(10,0))
        tk.Button(param_frame, text="Browse", 
                 command=lambda: self.browse_file(self.go_entry, 
                                                [("GAF Files", "*.gaf")])).grid(row=2, column=2, padx=5, pady=(10,0))
        
        tk.Label(param_frame, text=f"GO OBO File: {self.GO_OBO_PATH}").grid(row=3, column=0, sticky=tk.W, padx=5, pady=(10,0))
        
        tk.Label(param_frame, text="Essential Genes File (Excel):").grid(row=4, column=0, sticky=tk.W, padx=5, pady=(10,0))
        self.ess_entry = tk.Entry(param_frame, width=30)
        self.ess_entry.grid(row=4, column=1, padx=5, pady=(10,0))
        tk.Button(param_frame, text="Browse", 
                 command=lambda: self.browse_file(self.ess_entry, 
                                                [("Excel Files", "*.xlsx")])).grid(row=4, column=2, padx=5, pady=(10,0))

    def setup_search_button(self):
        # ----- Search Button ----- #
        search_button = tk.Button(self.root, text="Search", command=self.run_search)
        search_button.grid(row=2, column=0, pady=10)

    def setup_results_frame(self):
        # ----- Results Table with Scrollbars ----- #
        results_frame = tk.Frame(self.root)
        results_frame.grid(row=4, column=0, sticky="nsew", padx=10, pady=10)
        
        # Configure the frame to expand properly
        results_frame.rowconfigure(0, weight=1)
        results_frame.columnconfigure(0, weight=1)
        
        # Create a frame to hold both the treeview and scrollbars
        tree_frame = tk.Frame(results_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
    # Create Treeview with selectmode="extended" to allow multiple selections
        self.tree = ttk.Treeview(tree_frame, columns=(
            "Gene ID", "Protein", "Weight/kDa", "Length/aa", "Position", 
            "Matched Sequence", "P1' Residue", "Context", "Essential", "Located In", "GO Annotations"), 
            show="headings", selectmode="extended")
        
        
        # Set up column headings
        self.tree.heading("Gene ID", text="Gene ID", 
                        command=lambda: self.sort_treeview("Gene ID", False))
        self.tree.column("Gene ID", width=50, anchor=tk.W, stretch=False)
        
        self.tree.heading("Protein", text="Protein", 
                        command=lambda: self.sort_treeview("Protein", False))
        self.tree.column("Protein", width=300, anchor=tk.W, stretch=True)
        
        self.tree.heading("Weight/kDa", text="Weight/kDa", 
                        command=lambda: self.sort_treeview("Weight/kDa", False))
        self.tree.column("Weight/kDa", width=100, anchor=tk.CENTER, stretch=False)
        
        self.tree.heading("Length/aa", text="Length/aa", 
                        command=lambda: self.sort_treeview("Length/aa", False))
        self.tree.column("Length/aa", width=100, anchor=tk.CENTER, stretch=False)
        
        self.tree.heading("Position", text="Position", 
                        command=lambda: self.sort_treeview("Position", False))
        self.tree.column("Position", width=80, anchor=tk.CENTER, stretch=False)
        
        self.tree.heading("Matched Sequence", text="Matched Sequence", 
                        command=lambda: self.sort_treeview("Matched Sequence", False))
        self.tree.column("Matched Sequence", width=120, anchor=tk.W, stretch=False)
        
        self.tree.heading("P1' Residue", text="P1' Residue", 
                        command=lambda: self.sort_treeview("P1' Residue", False))
        self.tree.column("P1' Residue", width=80, anchor=tk.CENTER, stretch=False)
        
        self.tree.heading("Context", text="Context", 
                        command=lambda: self.sort_treeview("Context", False))
        self.tree.column("Context", width=200, anchor=tk.W, stretch=True)
        
        self.tree.heading("Essential", text="Essential", 
                        command=lambda: self.sort_treeview("Essential", False))
        self.tree.column("Essential", width=50, anchor=tk.CENTER, stretch=False)
        
        self.tree.heading("Located In", text="Located In", 
                        command=lambda: self.sort_treeview("Located In", False))
        self.tree.column("Located In", width=200, anchor=tk.W, stretch=True)
        
        self.tree.heading("GO Annotations", text="GO Annotations", 
                        command=lambda: self.sort_treeview("GO Annotations", False))
        self.tree.column("GO Annotations", width=200, anchor=tk.W, stretch=True)

        # Create vertical scrollbar
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        
        # Create horizontal scrollbar
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(xscrollcommand=hsb.set)
        
        # Grid layout for proper scrollbar functionality
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        # Configure grid weights
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)
    
        # Add binding for selection changes
        self.tree.bind("<<TreeviewSelect>>", self.update_selected_count)

    def setup_bottom_frame(self):
        # ----- Bottom Frame: Summary and Export ----- #
        bottom_frame = tk.Frame(self.root)
        bottom_frame.grid(row=5, column=0, pady=5)  # Changed to row 5
        
         # Add a label for selected items count
        self.selected_items_label = tk.Label(bottom_frame, text="Selected: 0", font=("Arial", 12, "bold"))
        self.selected_items_label.pack(side=tk.LEFT, padx=10)


        self.total_hits_label = tk.Label(bottom_frame, text="Total Hits Found: 0", font=("Arial", 12, "bold"))
        self.total_hits_label.pack(side=tk.LEFT, padx=10)
        
        # Add Structure Analysis button
        structure_button = tk.Button(bottom_frame, text="Structure Analysis", command=self.show_structure_analysis)
        structure_button.pack(side=tk.LEFT, padx=10)
        
        export_button = tk.Button(bottom_frame, text="Export to Excel", command=self.export_results)
        export_button.pack(side=tk.LEFT, padx=10)

    # Create a new frame for GO filtering
    def setup_go_filter_frame(self):
        go_filter_frame = tk.LabelFrame(self.root, text="GO Annotation Filters")
        go_filter_frame.grid(row=3, column=0, pady=5, padx=10, sticky="ew")  # Changed to row 3, full width
        
        # Create a grid layout within the frame
        go_filter_frame.columnconfigure(1, weight=1)  # Make the entry expandable
        
        # GO Filter entry
        tk.Label(go_filter_frame, text="GO Filter Expression:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.go_filter_entry = tk.Entry(go_filter_frame, width=40)
        self.go_filter_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        
        # Help text
        help_text = "Examples: membrane (Find all with 'membrane' in description) | GO:0005886 (Find exact GO ID) | membrane AND cytoplasm (Require both terms) | membrane OR transport (Either term) | membrane AND NOT nucleus (Exclude term)"
        help_label = tk.Label(go_filter_frame, text=help_text, justify=tk.LEFT, 
                            font=("Arial", 9), wraplength=800)  # Increased wraplength
        help_label.grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5, pady=5)
        
        # Filter controls in a separate frame
        filter_btn_frame = tk.Frame(go_filter_frame)
        filter_btn_frame.grid(row=0, column=2, padx=5, pady=5, sticky="e")
        
        # Add Filter button
        add_filter_button = tk.Button(filter_btn_frame, text="Apply Filter", command=self.apply_go_filter)
        add_filter_button.pack(side=tk.LEFT, padx=5)
        
        # Clear Filter button
        clear_filter_button = tk.Button(filter_btn_frame, text="Clear Filter", command=self.clear_go_filter)
        clear_filter_button.pack(side=tk.LEFT, padx=5)
        
        # Radio buttons in their own frame
        radio_frame = tk.Frame(go_filter_frame)
        radio_frame.grid(row=0, column=3, padx=5, pady=5, sticky="e")
        
        # Radio buttons for filter mode
        self.filter_mode = tk.StringVar(value="display")
        tk.Radiobutton(radio_frame, text="Filter Display Only", 
                    variable=self.filter_mode, value="display").pack(side=tk.LEFT, padx=5)
        tk.Radiobutton(radio_frame, text="Include in Search", 
                    variable=self.filter_mode, value="search").pack(side=tk.LEFT, padx=5)
   
   
    # ----- Helper Functions ----- #
    def get_default_genomes(self):
        if not os.path.exists(self.GENOME_FOLDER_PATH):
            os.makedirs(self.GENOME_FOLDER_PATH)
        return [f for f in os.listdir(self.GENOME_FOLDER_PATH) 
                if f.endswith((".fasta", ".fa", ".faa", ".txt"))]

    def update_fasta_entry(self, event):
        selected = self.genome_dropdown.get()
        if selected:
            self.fasta_file_entry.delete(0, tk.END)
            self.fasta_file_entry.insert(0, os.path.join(self.GENOME_FOLDER_PATH, selected))
            self.update_related_files()

    def browse_file(self, entry_widget, filetypes):
        file_path = filedialog.askopenfilename(filetypes=filetypes)
        if file_path:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, file_path)
            self.genome_dropdown.set("")
            self.update_related_files()

    def update_related_files(self):
        """Automatically update GO and Essential file entries based on the genome filename."""
        fasta_path = self.fasta_file_entry.get().lower()
        if "ecoli" in fasta_path:
            self.go_entry.delete(0, tk.END)
            self.go_entry.insert(0, os.path.join("go_files", "ecocyc.gaf"))
            self.ess_entry.delete(0, tk.END)
            self.ess_entry.insert(0, os.path.join("go_files", "essential_ecoli.xlsx"))
        elif "human" in fasta_path:
            self.go_entry.delete(0, tk.END)
            self.go_entry.insert(0, os.path.join("go_files", "goa_human.gaf"))
            self.ess_entry.delete(0, tk.END)
            self.ess_entry.insert(0, "")

    def sort_treeview(self, col, reverse):
        """Sort Treeview contents when header is clicked."""
        l = [(self.tree.set(k, col), k) for k in self.tree.get_children('')]
        try:
            l.sort(key=lambda t: float(t[0]), reverse=reverse)
        except ValueError:
            l.sort(reverse=reverse)
        for index, (val, k) in enumerate(l):
            self.tree.move(k, '', index)
        self.tree.heading(col, command=lambda: self.sort_treeview(col, not reverse))
    
    #method to update the selection count and bind it to the Treeview:
    def update_selected_count(self, event=None):
        """Update the count of selected items in the tree"""
        selected_count = len(self.tree.selection())
        self.selected_items_label.config(text=f"Selected: {selected_count}")    

    # ----- Main Functionality ----- #
    def run_search(self):
        fasta_path = self.fasta_file_entry.get()
        if not fasta_path or not os.path.exists(fasta_path):
            messagebox.showerror("Error", "Please provide a valid genome file.")
            return
        regex_pat = self.regex_entry.get().strip()
        if not regex_pat:
            messagebox.showerror("Error", "Please enter a regex pattern.")
            return
        p1_input = self.p1_entry.get().strip()
        p1_regex = f"[{p1_input.replace(',', '')}]" if p1_input else ""
        at_cterm = self.at_cterm_var.get()
        close_cterm = self.close_cterm_var.get()
        
        try:
            hits = self.sequence_analyzer.scan_fasta_for_regex(fasta_path, regex_pat, p1_regex, at_cterm, close_cterm)
            if not hits:
                messagebox.showinfo("Results", "No hits found.")
                self.current_hits = []
                self.filtered_hits = []  # Also clear filtered hits
                self.update_table([])
                self.total_hits_label.config(text="Total Hits Found: 0")
                return
            
            go_path = self.go_entry.get().strip()
            obo_path = self.GO_OBO_PATH
            if go_path and os.path.exists(go_path) and os.path.exists(obo_path):
                go_ann = self.sequence_analyzer.load_go_annotations(go_path)
                go_trans = self.sequence_analyzer.load_go_translations(obo_path)
                hits = self.sequence_analyzer.annotate_hits(hits, go_ann, go_trans)
            else:
                for hit in hits:
                    hit["go_terms"] = {}

            # Apply GO filter if one is set for search
            if hasattr(self, 'current_go_filter') and self.current_go_filter and self.filter_mode.get() == "search":
                hits = self.go_filter.filter_hits_by_expression(hits, self.current_go_filter)
                if not hits:
                    messagebox.showinfo("Results", "No hits found after applying GO filter.")
                    self.current_hits = []
                    self.filtered_hits = []  # Also clear filtered hits
                    self.update_table([])
                    self.total_hits_label.config(text="Total Hits Found: 0")
                    return        
            
            ess_path = self.ess_entry.get().strip()
            if "ecoli" in fasta_path.lower() and ess_path and os.path.exists(ess_path):
                hits = self.sequence_analyzer.mark_essential_genes(hits, ess_path)
            else:
                for hit in hits:
                    hit["Essential"] = "N/A"
            
            hits = self.sequence_analyzer.get_reduced_hits(hits)
            hits.sort(key=lambda x: x["gene_id"])
            self.current_hits = hits
            self.filtered_hits = hits.copy()  # Initialize filtered hits as a copy of all hits
            self.update_table(self.filtered_hits)
            self.total_hits_label.config(text=f"Total Hits Found: {len(hits)}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def update_table(self, data):
        """Update the table with data and maintain a mapping of tree item IDs to hit indices."""
        # Clear the existing mapping and tree
        self.tree_id_to_hit_index = {}
        for row in self.tree.get_children():
            self.tree.delete(row)
            
        for i, hit in enumerate(data):
            go_text = ""
            if hit.get("go_terms"):
                parts = []
                for qualifier, terms in hit["go_terms"].items():
                    for go_id, term in terms:
                        parts.append(f"{term}")
                go_text = "; ".join(parts)
            
            # Extract located_in annotations
            located_in_text = "N/A"
            if hit.get("go_terms") and "located_in" in hit["go_terms"]:
                located_in_parts = []
                for go_id, term in hit["go_terms"]["located_in"]:
                    located_in_parts.append(term)
                located_in_text = "; ".join(located_in_parts)
            
            item_id = self.tree.insert("", "end", values=(
                hit["gene_id"],
                hit["protein"],
                f"{hit['weight']/1000:.1f}",
                hit["protein_size"],
                hit["occurrences"][0][0]+1,
                hit["sequence"][0],
                hit["p1prime"][0],
                hit["match_snippets"][0],
                hit.get("Essential", "N/A"),
                located_in_text,
                go_text
            ))
            
            # Store the mapping from tree item ID to hit index
            self.tree_id_to_hit_index[item_id] = i

    def export_results(self):
        if not self.current_hits:
            messagebox.showerror("Error", "No results to export.")
            return
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", 
                                                filetypes=[("Excel Files", "*.xlsx")])
        if not file_path:
            return
        try:
            success, message = self.sequence_analyzer.export_hits_to_excel(self.current_hits, file_path)
            messagebox.showinfo("Success", message)
        except Exception as e:
            messagebox.showerror("Error", str(e))

   

    def show_structure_analysis(self):
        """Open a new window with structure analysis for selected items"""
        selected_items = self.tree.selection()
        
        if not selected_items:
            messagebox.showwarning("Selection Required", "Please select at least one item for structure analysis.")
            return
        
        # Get the selected hits from filtered_hits based on tree selection and our mapping
        selected_hits = []
        for item_id in selected_items:
            if item_id in self.tree_id_to_hit_index:
                hit_index = self.tree_id_to_hit_index[item_id]
                if hit_index < len(self.filtered_hits):  # Safety check
                    selected_hits.append(self.filtered_hits[hit_index])
        
        # Create a loading indicator window
        loading_window = tk.Toplevel(self.root)
        loading_window.title("Processing")
        loading_window.geometry("300x150")
        loading_window.grab_set()  # Make the window modal
        
        # Main message label
        loading_label = tk.Label(loading_window, 
                            text="Calculating structure analysis...\nThis may take some time.", 
                            pady=20)
        loading_label.pack()
        
        # Counter label
        counter_label = tk.Label(loading_window, 
                                text=f"Processing: 0/{len(selected_hits)}", 
                                font=("Arial", 10, "bold"))
        counter_label.pack()
        
        # Define a progress callback function
        def update_progress(current_index):
            counter_label.config(text=f"Processing: {current_index+1}/{len(selected_hits)}")
            loading_window.update()
        
        # Update the loading window to show it
        loading_window.update()
        
        # Determine which organism to use based on the selected genome file
        genome_path = self.fasta_file_entry.get().lower()
        organism = "human" if "human" in genome_path else "ecoli"
        
        # Use the structure analyzer to analyze the hits with progress callback
        analyzed_hits = self.structure_analyzer.analyze_proteins(selected_hits, organism, progress_callback=update_progress)
        
        # Close the loading window
        loading_window.destroy()

        # Store selected hits for later use
        self.current_structure_hits = analyzed_hits
        
        # Create a new window for structure analysis
        self.structure_window = tk.Toplevel(self.root)
        self.structure_window.title("Structure Analysis")
        self.structure_window.geometry("1200x800")
        self.structure_window.grab_set()  # Make the window modal
        
        # Set up the structure analysis table
        self.setup_structure_table(self.structure_window, analyzed_hits)

    def setup_structure_table(self, window, selected_hits):
        """Create the structure analysis table in the new window"""
        # Add title label
        title_label = tk.Label(window, text=f"Structure Analysis for {len(selected_hits)} Selected Items", 
                            font=("Arial", 12, "bold"))
        title_label.pack(side=tk.TOP, pady=10)
        
        # Add secondary structure legend
        legend_text = "Secondary Structure Legend: H = Alpha helix (4-12), B = Isolated beta-bridge residue, E = Strand, G = 3-10 helix, I = Pi helix, T = Turn, S = Bend, - = None"
        legend_label = tk.Label(window, text=legend_text, font=("Arial", 11), wraplength=980)
        legend_label.pack(side=tk.TOP, pady=5)

        # Add secondary structure legend
        legend_text = "Double Click Link to Uniprot"
        legend_label = tk.Label(window, text=legend_text, font=("Arial", 11), wraplength=980)
        legend_label.pack(side=tk.TOP, pady=5)
        
        # Create a frame for the table
        frame = tk.Frame(window)
        frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create the treeview with structure analysis columns
        columns = ("Gene ID",  "Protein", "Uniprot Link", "Essential", "Total SASA", 
                "Secondary Structure", "Residues within 5Å", 
                "Residues within 10Å", "Residues within 15Å")
        
        tree = ttk.Treeview(frame, columns=columns, show="headings")
        
        # Define column headings and properties
        tree.heading("Gene ID", text="Gene ID", command=lambda: self.sort_window_treeview(tree, "Gene ID", False))
        tree.column("Gene ID", width=50, anchor=tk.W)
        
        tree.heading("Protein", text="Protein", command=lambda: self.sort_window_treeview(tree, "Protein", False))
        tree.column("Protein", width=300, anchor=tk.W)

        tree.heading("Uniprot Link", text="Uniprot Link")
        tree.column("Uniprot Link", width=100, anchor=tk.CENTER)
        
        tree.heading("Essential", text="Essential", command=lambda: self.sort_window_treeview(tree, "Essential", False))
        tree.column("Essential", width=50, anchor=tk.CENTER)
        
        tree.heading("Total SASA", text="Total SASA", command=lambda: self.sort_window_treeview(tree, "Total SASA", False))
        tree.column("Total SASA", width=100, anchor=tk.CENTER)
        
        tree.heading("Secondary Structure", text="Secondary Structure", 
                    command=lambda: self.sort_window_treeview(tree, "Secondary Structure", False))
        tree.column("Secondary Structure", width=100, anchor=tk.W)
        
        tree.heading("Residues within 5Å", text="Residues within 5Å", 
                    command=lambda: self.sort_window_treeview(tree, "Residues within 5Å", False))
        tree.column("Residues within 5Å", width=100, anchor=tk.CENTER)
        
        tree.heading("Residues within 10Å", text="Residues within 10Å", 
                    command=lambda: self.sort_window_treeview(tree, "Residues within 10Å", False))
        tree.column("Residues within 10Å", width=100, anchor=tk.CENTER)
        
        tree.heading("Residues within 15Å", text="Residues within 15Å", 
                    command=lambda: self.sort_window_treeview(tree, "Residues within 15Å", False))
        tree.column("Residues within 15Å", width=100, anchor=tk.CENTER)
        
        # Add scrollbars
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vsb.set)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(xscrollcommand=hsb.set)
        
        # Place elements in the grid
        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        # Configure the grid
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        
        # Add double-click binding to open Uniprot links
        tree.bind("<Double-1>", lambda event: self.open_uniprot_link(event, tree))
        
        # Fill the table with structure data
        for hit in selected_hits:
            # Format secondary structure to a continuous string without commas
            sec_struct = self.structure_analyzer.format_secondary_structure(hit.get("secondary_structure", {}))
            
            # Get Uniprot link text - always show link if we have a valid uniprot_id
            uniprot_id = hit.get("uniprot_id")
            # Show "N/A" if we couldn't find the Uniprot ID at all
            if uniprot_id and uniprot_id != "Not found":
                # Remove any "sp|" prefix if present
                uniprot_id = uniprot_id.split('|')[-1]
                uniprot_link = f"Uniprot {uniprot_id}"
            else:
                uniprot_link = "N/A"
            
            # Insert the row with all values
            tree.insert("", "end", values=(
                hit["gene_id"],
                hit["protein"],
                uniprot_link,
                hit.get("Essential", "N/A"),
                f"{hit.get('total_sasa', 'N/A'):.2f}" if isinstance(hit.get('total_sasa'), (float, int)) else "N/A",
                sec_struct,
                hit.get("residues_within_5.0A", "N/A"),
                hit.get("residues_within_10.0A", "N/A"),
                hit.get("residues_within_15.0A", "N/A")
            ))
        
        # Add export button
        export_frame = tk.Frame(window)
        export_frame.pack(pady=10)
        
        export_button = tk.Button(export_frame, text="Export Structure Data", 
                                command=lambda: self.export_structure_data(selected_hits))
        export_button.pack()

    def open_uniprot_link(self, event, tree):
        """Open Uniprot link when double-clicked"""
        # Get the clicked item
        item = tree.identify('item', event.x, event.y)
        if not item:
            return
        
        # Get the column that was clicked
        column = tree.identify_column(event.x)
        column_index = int(column.replace('#', '')) - 1
        column_name = tree.column(column, "id")
        
        # Check if the Uniprot Link column was clicked
        if column_name != "Uniprot Link":
            return
        
        # Get the item values
        row_values = tree.item(item, 'values')
        if not row_values:
            return
            
        # Get the gene_id from the first column
        gene_id = row_values[0]
        
        # Find the corresponding hit from current_structure_hits
        for hit in self.current_structure_hits:
            if hit["gene_id"] == gene_id:
                uniprot_id = hit.get("uniprot_id")
                if uniprot_id and uniprot_id != "Not found":
                    # Remove any "sp|" prefix if present
                    uniprot_id = uniprot_id.split('|')[-1]
                    url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}/entry"
                    webbrowser.open_new_tab(url)
                break

    def apply_go_filter(self):
        """Apply GO filter to current hits or include in search"""
        filter_expression = self.go_filter_entry.get().strip()
        # Remove double quotes using regular expressions
        filter_expression = re.sub(r'"', '', filter_expression)

        # Replace commas with ' OR '
        filter_expression = filter_expression.replace(',', ' OR ')

        if not filter_expression:
            return
        
        # If no hits yet and we're in display mode, show a message
        if not self.current_hits and self.filter_mode.get() == "display":
            messagebox.showinfo("No Data", "Please run a search first before applying display filters.")
            return
            
        if self.filter_mode.get() == "display" and self.current_hits:
            # Filter the current display
            self.filtered_hits = self.go_filter.filter_hits_by_expression(self.current_hits, filter_expression)
            self.update_table(self.filtered_hits)
            self.total_hits_label.config(text=f"Filtered Hits: {len(self.filtered_hits)} of {len(self.current_hits)}")
        else:
            # Store the filter for use during search
            self.current_go_filter = filter_expression
            messagebox.showinfo("Filter Set", f"GO filter '{filter_expression}' will be applied in the next search.")

    def sort_window_treeview(self, tree, col, reverse):
        """Sort Treeview contents in the structure window"""
        l = [(tree.set(k, col), k) for k in tree.get_children('')]
        try:
            l.sort(key=lambda t: float(t[0]), reverse=reverse)
        except ValueError:
            l.sort(reverse=reverse)
        for index, (val, k) in enumerate(l):
            tree.move(k, '', index)
        tree.heading(col, command=lambda: self.sort_window_treeview(tree, col, not reverse))
    
    def export_structure_data(self, selected_hits):
        """Export the structure analysis data to Excel"""
        if not selected_hits:
            messagebox.showinfo("Info", "No data to export.")
            return
        
        file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", 
                                            filetypes=[("Excel Files", "*.xlsx")])
        if not file_path:
            return
        
        success, message = self.sequence_analyzer.export_structure_data_to_excel(selected_hits, file_path)
        if success:
            messagebox.showinfo("Success", message)
        else:
            messagebox.showerror("Error", message)

        
    def clear_go_filter(self):
        """Clear any applied GO filters"""
        self.go_filter_entry.delete(0, tk.END)
        self.current_go_filter = ""
        
        # If we have hits, restore the full display
        if self.current_hits:
            self.filtered_hits = self.current_hits.copy()  # Reset filtered hits to all hits
            self.update_table(self.current_hits)
            self.total_hits_label.config(text=f"Total Hits Found: {len(self.current_hits)}")

# Entry point for the UI
if __name__ == "__main__":
    root = tk.Tk()
    app = SequenceAnalyzerUI(root)
    root.mainloop()