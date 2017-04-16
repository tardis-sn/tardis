(function(){
  'use strict';
  /**
   * Tag input directive for use in list fields
   * based on http://codepen.io/webmatze/pen/isuHh
   */
  angular
      .module('tardisApp')
      .directive('tagInput', tagInput);

  function tagInput() {
    var directive = {
      restrict: 'E',
      templateUrl: '/static/app/input-data/tag.tmpl.html',
      scope: {
        inputTags: '=taglist',        
      },    
      link: link,
      controller: tagInputCtrl,
      controllerAs: 'vm',
      bindToController: true
    }
    return directive;

    function link (scope, elem, attrs){
      scope.vm.defaultWidth = 200;
      scope.vm.tagText = '';
      scope.vm.placeholder = attrs.placeholder;
      scope.vm.tagArray = tagArray;
      scope.vm.addTag = addTag;
      scope.vm.deleteTag = deleteTag;

      scope.$watch('tagText', watch);
      elem.bind('keydown', keydown);
      return elem.bind('keyup', keyup);

      function tagArray() {
        if (scope.vm.inputTags === undefined) {
          return [];
        }
        return scope.vm.inputTags.split(',').filter(function(tag) {
          return tag !== "";
        });
      }
          
      function addTag () {
        var tagArray;
        if (scope.vm.tagText.length === 0) {
          console.log('lalalla');
          return;
        }
        tagArray = scope.vm.tagArray();
        //evaluate if new is already in array
        var inx = tagArray.indexOf(scope.vm.tagText);
        if (inx == -1){
          tagArray.push(scope.vm.tagText);
          scope.vm.inputTags = tagArray.join(',');
          return scope.vm.tagText = "";
        }else{
          window.alert('repeated value');
        }
      }
      
      function deleteTag(key) {
        var tagArray;
        tagArray = scope.vm.tagArray();
        if (tagArray.length > 0 && scope.vm.tagText.length === 0 && key === undefined) {
          tagArray.pop();
        } else {
          if (key !== undefined) {
            tagArray.splice(key, 1);
          }
        }
        return scope.vm.inputTags = tagArray.join(',');
      }
       
      function watch(newVal, oldVal) {
        var tempEl;
        if (!(newVal === oldVal && newVal === undefined)) {
          tempEl = $("<span>" + newVal + "</span>").appendTo("body");
          scope.vm.inputWidth = tempEl.width() + 5;
          if (scope.vm.inputWidth < scope.vm.defaultWidth) {
            scope.vm.inputWidth = scope.vm.defaultWidth;
          }
          return tempEl.remove();
        }
      }
      
      function keydown(e) {
        var key;
        key = e.which;
        if (key === 9 || key === 13) {
          e.preventDefault();
        }
        if (key === 8 && scope.vm.tagText.length === 0) {
          return scope.$apply(scope.vm.deleteTag);
        }
      }
        
      function keyup(e) {
        var key;
        key = e.which;
        if (key === 9 || key === 13 || key === 188) {
          e.preventDefault();
          return scope.$apply(scope.vm.addTag);
        }
      }

    }

  }

  tagInputCtrl.$inject = ['$scope'];

  function tagInputCtrl($scope){
    var vm = this;
  }
})();