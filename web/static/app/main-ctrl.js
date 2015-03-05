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
  vm.model = {};
  vm.model.struture = 'file';
  vm.model.abundance = 'uniform';
  vm.model.density = 'branch85_w7';

  vm.changeTab = changeTab;
  vm.nextTab = nextTab;
  vm.prevTab = prevTab;
  vm.nextButtonName = 'next'

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