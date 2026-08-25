use chem_core::molecule::Molecule;
use egui::Ui;

pub fn show_molecule_info(ui: &mut Ui, mol: &Molecule, smiles: &str, name: &str) {
    ui.group(|ui| {
        ui.heading("Molecule Information");

        ui.separator();

        egui::Grid::new("mol_info_grid")
            .num_columns(2)
            .spacing([10.0, 5.0])
            .striped(true)
            .show(ui, |ui| {
                ui.label(egui::RichText::new("Name:").strong());
                ui.label(name);
                ui.end_row();

                ui.label(egui::RichText::new("SMILES:").strong());
                ui.label(egui::RichText::new(smiles).code());
                ui.end_row();

                ui.label(egui::RichText::new("Formula:").strong());
                ui.label(mol.formula());
                ui.end_row();

                ui.label(egui::RichText::new("Molecular Weight:").strong());
                ui.label(format!("{:.2} g/mol", mol.molecular_weight()));
                ui.end_row();

                ui.label(egui::RichText::new("Atoms:").strong());
                ui.label(format!("{}", mol.num_atoms()));
                ui.end_row();

                ui.label(egui::RichText::new("Bonds:").strong());
                ui.label(format!("{}", mol.num_bonds()));
                ui.end_row();
            });
    });
}

pub fn show_atom_list(ui: &mut Ui, mol: &Molecule) {
    ui.group(|ui| {
        ui.label(egui::RichText::new("Atoms").strong().size(14.0));

        egui::ScrollArea::vertical()
            .id_salt("atom_scroll")
            .max_height(200.0)
            .show(ui, |ui| {
                egui::Grid::new(format!("atom_list_{:p}", mol as *const _))
                    .num_columns(4)
                    .spacing([10.0, 5.0])
                    .striped(true)
                    .show(ui, |ui| {
                        ui.label(egui::RichText::new("ID").strong());
                        ui.label(egui::RichText::new("Element").strong());
                        ui.label(egui::RichText::new("Degree").strong());
                        ui.label(egui::RichText::new("Aromatic").strong());
                        ui.end_row();

                        for (i, atom) in mol.atoms().iter().enumerate() {
                            ui.label(format!("{}", i));
                            ui.label(atom.element().symbol());
                            ui.label(format!("{}", mol.degree(i)));
                            ui.label(if atom.is_aromatic() { "Yes" } else { "No" });
                            ui.end_row();
                        }
                    });
            });
    });
}

pub fn show_bond_list(ui: &mut Ui, mol: &Molecule) {
    ui.group(|ui| {
        ui.label(egui::RichText::new("Bonds").strong().size(14.0));

        egui::ScrollArea::vertical()
            .id_salt("bond_scroll")
            .max_height(200.0)
            .show(ui, |ui| {
                egui::Grid::new(format!("bond_list_{:p}", mol as *const _))
                    .num_columns(4)
                    .spacing([10.0, 5.0])
                    .striped(true)
                    .show(ui, |ui| {
                        ui.label(egui::RichText::new("ID").strong());
                        ui.label(egui::RichText::new("Atoms").strong());
                        ui.label(egui::RichText::new("Order").strong());
                        ui.label(egui::RichText::new("Aromatic").strong());
                        ui.end_row();

                        for (i, bond) in mol.bonds().iter().enumerate() {
                            ui.label(format!("{}", i));
                            ui.label(format!("{} - {}", bond.atom1(), bond.atom2()));
                            ui.label(format!("{:?}", bond.order()));
                            ui.label(if bond.is_aromatic() { "Yes" } else { "No" });
                            ui.end_row();
                        }
                    });
            });
    });
}

pub fn molecule_compact(ui: &mut Ui, smiles: &str, name: &str, formula: &str) {
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(name).strong());
        ui.separator();
        ui.label(egui::RichText::new(smiles).code().small());
        ui.separator();
        ui.label(formula);
    });
}
