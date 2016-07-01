#include "dialog_modelvars.h"
#include "ui_dialog_modelvars.h"
#include <QDoubleValidator>

dialog_modelvars::dialog_modelvars(QWidget *parent) 
	: QDialog(parent),
	ui(new Ui::dialog_modelvars_class)
{
	ui->setupUi(this);
	//Initialize boxvalidator
	ui->lineStackSpacing->setValidator(new QDoubleValidator(0, 1000, 4, this));
	ui->lineDiameter->setValidator(new QDoubleValidator(0, 1000, 4, this));
}

dialog_modelvars::~dialog_modelvars()
{
	delete ui;
}

//Custom functions
void dialog_modelvars::on_okButton_clicked()
{
	dialogOK = true;
	m_nr_atoms = int(ui->spbNrAtoms->value());
	m_stackspacing = double(ui->lineStackSpacing->text().toDouble());
	m_diameter = double(ui->lineDiameter->text().toDouble());
	dialog_modelvars::accept();
}

void dialog_modelvars::on_cancelButton_clicked()
{
	dialogOK = false;
	dialog_modelvars::reject();
}
