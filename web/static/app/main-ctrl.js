'use strict';

angular
  .module('tardisApp')
  .controller('MainCtrl', MainCtrl);

function MainCtrl (){
  var vm = this;

  vm.title = 'pruba de esto'
  vm.tabs = [
    {
      'name': 'Supernova',
      'template': 'input-data/supernova-partial.html'
    },
    {
      'name': 'Atom data & Spectrum',
      'template': 'input-data/atom-spectrum-partial.html'
    },
    {
      'name': 'Plasma',
      'template': 'input-data/plasma-partial.html'
    },
    {
      'name': 'montecarlo',
      'template': 'input-data/montecarlo-partial.html'
    },
    {
      'name': 'Model',
      'template': 'input-data/model-partial.html'
    }
  ]
  vm.tabIndex = 0;
  vm.CurrentTab = vm.tabs[vm.tabIndex];
  // define the ng-models object for input data
  vm.data = {};
  vm.data.supernova = {};
  vm.data.atom_data = {};
  vm.data.spectrum = {};
  vm.data.plasma = {};
  vm.data.plasma.nlte = {};
  vm.data.montecarlo = {};
  vm.data.montecarlo.black_body_sampling = {};
  vm.data.model = {};
  vm.data.model.struture = {};
  vm.data.model.struture.type = 'file';
  vm.data.model.struture.velocity = {}
  vm.data.model.struture.density = {};
  vm.data.model.struture.density.type = 'branch85_w7';
  vm.data.model.abundance = {};
  vm.data.model.abundance.type = 'uniform';

  vm.changeTab = changeTab;
  vm.nextTab = nextTab;
  vm.prevTab = prevTab;

  function changeTab(index, tab) {
    vm.tabIndex = index;
    vm.CurrentTab = tab;
  };

  function nextTab(){
    if (vm.tabIndex < vm.tabs.length - 1) {
      vm.tabIndex++;
      vm.CurrentTab = vm.tabs[vm.tabIndex];
    }
  };

  function prevTab(){
    if (vm.tabIndex > 0) {
      vm.tabIndex--;
      vm.CurrentTab = vm.tabs[vm.tabIndex];
    }
  };

}